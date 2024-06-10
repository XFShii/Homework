import glob
import json
import logging
import os
import re
import dpdata
import numpy as np
from monty.serialization import dumpfn, loadfn
from pymatgen.core.structure import Structure
from pymatgen.core.surface import SlabGenerator

from apex.core.calculator.lib import abacus_utils
from apex.core.calculator.lib import vasp_utils
from apex.core.property.Property import Property
from apex.core.refine import make_refine
from apex.core.reproduce import make_repro, post_repro
from dflow.python import upload_packages
upload_packages.append(__file__)

class DecohesionEnergy(Property):
    def __init__(self, parameter, inter_param=None):
        '''
               core parameter for make_confs[POSCAR]：
               min_slab_size  max_vacuum_size  vacuum_size_step
               miller_index
        '''
        parameter["reproduce"] = parameter.get("reproduce", False)
        self.reprod = parameter["reproduce"]
        if not self.reprod:
            if not ("init_from_suffix" in parameter and "output_suffix" in parameter):
                self.min_slab_size = parameter["min_slab_size"]
                parameter["pert_xz"] = parameter.get("pert_xz", 0.01)
                self.pert_xz = parameter["pert_xz"]
                parameter["max_vacuum_size"] = parameter.get("max_vacuum_size", 15)
                self.max_vacuum_size = parameter["max_vacuum_size"]
                parameter["vacuum_size_step"]=parameter.get("vacuum_size_step", 1)
                self.vacuum_size_step = parameter["vacuum_size_step"]
                self.miller_index = tuple(parameter["miller_index"])
            parameter["cal_type"] = parameter.get("cal_type", "relaxation")
            default_cal_setting = {
                "relax_pos": False,
                "relax_shape": False,
                "relax_vol": False,
            }
        else:
            parameter["cal_type"] = "static"
            self.cal_type = parameter["cal_type"]
            default_cal_setting = {
                "relax_pos": False,
                "relax_shape": False,
                "relax_vol": False,
            }
            parameter["init_from_suffix"] = parameter.get("init_from_suffix", "00")
            self.init_from_suffix = parameter["init_from_suffix"]
        self.cal_type = parameter["cal_type"]
        parameter["cal_setting"] = parameter.get("cal_setting", default_cal_setting)
        for key in default_cal_setting:
            parameter["cal_setting"].setdefault(key, default_cal_setting[key])
        self.cal_setting = parameter["cal_setting"]
        self.parameter = parameter
        self.inter_param = inter_param if inter_param != None else {"type": "vasp"}

    def make_confs(self, path_to_work, path_to_equi, refine=False):
        path_to_work = os.path.abspath(path_to_work)
        if os.path.exists(path_to_work):
            logging.warning("%s already exists" % path_to_work)
        else:
            os.makedirs(path_to_work)
        path_to_equi = os.path.abspath(path_to_equi)

        task_list = []
        cwd = os.getcwd()

        if self.reprod:
            print("surface reproduce starts")
            if "init_data_path" not in self.parameter:
                raise RuntimeError("please provide the initial data path to reproduce")
            init_data_path = os.path.abspath(self.parameter["init_data_path"])
            task_list = make_repro(
                self.inter_param,
                init_data_path,
                self.init_from_suffix,
                path_to_work,
                self.parameter.get("reprod_last_frame", True),
            )

        else:
            if refine:
                logging.info("surface refine starts")
                task_list = make_refine(
                    self.parameter["init_from_suffix"],
                    self.parameter["output_suffix"],
                    path_to_work,
                )
                # record miller
                init_from_path = re.sub(
                    self.parameter["output_suffix"][::-1],
                    self.parameter["init_from_suffix"][::-1],
                    path_to_work[::-1],
                    count=1,
                )[::-1]
                task_list_basename = list(map(os.path.basename, task_list))

                for ii in task_list_basename:
                    init_from_task = os.path.join(init_from_path, ii)
                    output_task = os.path.join(path_to_work, ii)
                    os.chdir(output_task)
                    if os.path.isfile("decohesion_energy.json"):
                        os.remove("decohesion_energy.json")
                    if os.path.islink("decohesion_energy.json"):
                        os.remove("decohesion_energy.json")
                    os.symlink(
                        os.path.relpath(os.path.join(init_from_task, "decohesion_energy.json")),
                        "decohesion_energy.json",)
            else:
                if self.inter_param["type"] == "abacus":
                    CONTCAR = abacus_utils.final_stru(path_to_equi)
                    POSCAR = "STRU"
                else:
                # refine = false && reproduce = false && self.inter_param["type"]== "vasp"
                    CONTCAR = "CONTCAR"
                    POSCAR = "POSCAR"

                equi_contcar = os.path.join(path_to_equi, CONTCAR)
                if not os.path.exists(equi_contcar):
                    raise RuntimeError("please do relaxation first")

                if self.inter_param["type"] == "abacus":
                    stru = dpdata.System(equi_contcar, fmt="stru")
                    stru.to("contcar", "CONTCAR.tmp")
                    ptypes = vasp_utils.get_poscar_types("CONTCAR.tmp")
                    ss = Structure.from_file("CONTCAR.tmp")
                    os.remove("CONTCAR.tmp")
                else:
                    ptypes = vasp_utils.get_poscar_types(equi_contcar)
                    # element type 读取 vasp 第五行
                    ss = Structure.from_file(equi_contcar)
                # gen POSCAR of Slab
                all_slabs = []
                vacuum = []
                num = 0
                while self.vacuum_size_step * num <= self.max_vacuum_size:
                    vacuum_size = self.vacuum_size_step * num
                    gen = SlabGenerator(
                        ss,
                        self.miller_index,
                        self.min_slab_size,
                        vacuum_size,
                    )
                    slabs = gen.get_slabs()
                    num = num + 1
                    all_slabs.extend(slabs)
                    vacuum.append(vacuum_size)

                os.chdir(path_to_work)
                if os.path.exists(POSCAR):
                    os.remove(
                        POSCAR)
                os.symlink(os.path.relpath(equi_contcar), POSCAR)
                for ii in range(len(all_slabs)):
                    output_task = os.path.join(path_to_work, "task.%06d" % ii)
                    os.makedirs(output_task, exist_ok=True)
                    os.chdir(output_task)
                    for jj in [
                        "INCAR",
                        "POTCAR",
                        "POSCAR",
                        "conf.lmp",
                        "in.lammps",
                        "STRU",
                    ]:
                        if os.path.exists(jj):
                            os.remove(jj)
                    task_list.append(output_task)
                    logging.info(
                        "# %03d generate " % ii,
                        output_task,
                        " \t %d atoms" % len(all_slabs[ii].sites),
                    )
                    all_slabs[ii].to("POSCAR.tmp", "POSCAR")
                    vasp_utils.regulate_poscar("POSCAR.tmp", "POSCAR")
                    vasp_utils.sort_poscar("POSCAR", "POSCAR", ptypes)
                    vasp_utils.perturb_xz("POSCAR", "POSCAR", self.pert_xz)
                    if self.inter_param["type"] == "abacus":
                        abacus_utils.poscar2stru("POSCAR", self.inter_param, "STRU")
                        os.remove("POSCAR")
                    decohesion_energy = {"miller_index": all_slabs[ii].miller_index, "vacuum_size":vacuum[ii]}
                    dumpfn(decohesion_energy, "decohesion_energy.json", indent=4)
        os.chdir(cwd)
        return task_list

    def post_process(self, task_list):
        pass

    def task_type(self):
        return self.parameter["type"]

    def task_param(self):
        return self.parameter

    def _compute_lower(self, output_file, all_tasks, all_res) -> [dict, str]:
        output_file = os.path.abspath(output_file)
        res_data = {}
        ptr_data = os.path.dirname(output_file) + "\n"
        if not self.reprod:
            equi_path = os.path.abspath(
                os.path.join(
                    os.path.dirname(output_file), "../relaxation/relax_task"
                )
            )
            equi_result = loadfn(os.path.join(equi_path, "result.json"))
            equi_epa = equi_result["energies"][-1] / np.sum(
                equi_result["atom_numbs"]
            )
            vacuum_size_step = loadfn(os.path.join(os.path.dirname(output_file), "param.json"))["vacuum_size_step"]
            ptr_data += ("Miller Index: " + str(loadfn(os.path.join(os.path.dirname(output_file), "param.json"))["miller_index"]) + "\n")
            ptr_data += "Vacuum_size(e-10 m):\tDecohesion_E(J/m^2) Decohesion_S(Pa)\n"
            pre_evac = 0
            for ii in all_tasks:
                task_result = loadfn(os.path.join(ii, "result_task.json"))
                natoms = np.sum(task_result["atom_numbs"])
                epa = task_result["energies"][-1] / natoms
                AA = np.linalg.norm(
                    np.cross(task_result["cells"][0][0], task_result["cells"][0][1])
                )

                structure_dir = os.path.basename(ii)
                Cf = 1.60217657e-16 / (1e-20 * 2) * 0.001
                evac = (task_result["energies"][-1] - equi_epa * natoms) / AA * Cf
                vacuum_size = loadfn(os.path.join(ii, "decohesion_energy.json"))["vacuum_size"]
                stress = (evac - pre_evac) / vacuum_size_step * 1e10

                ptr_data += "%-30s   % 7.3f     %10.3e\n" % (
                    str(vacuum_size) + "-" + structure_dir + ":",
                    evac,
                    stress,
                )
                res_data[str(vacuum_size) + "_" + structure_dir] = [
                    evac,
                    stress,
                    vacuum_size,
                ]
                pre_evac = evac

        else:
            if "init_data_path" not in self.parameter:
                raise RuntimeError("please provide the initial data path to reproduce")
            init_data_path = os.path.abspath(self.parameter["init_data_path"])
            res_data, ptr_data = post_repro(
                init_data_path,
                self.parameter["init_from_suffix"],
                all_tasks,
                ptr_data,
                self.parameter.get("reprod_last_frame", True),
            )

        with open(output_file, "w") as fp:
            json.dump(res_data, fp, indent=4)

        return res_data, ptr_data


