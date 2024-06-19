import plotly.graph_objs as go
import pandas as pd
from apex.core.lib.utils import round_format, round_2d_format
def plotly_graph(res_data: dict, name: str, **kwargs):
    decohesion_e = [values[0] for values in res_data.values()]
    stress = [values[1] for values in res_data.values()]
    vacuum_size = [values[2] for values in res_data.values()]
    vacuum_size = [str(item) for item in vacuum_size]
    df = pd.DataFrame({
        "separation distance (A)": vacuum_size,
        "Decohesion energy (J/m^2)": decohesion_e,
        "Decohesion stress (GPa)": [ s / 1e9 for s in stress],
    })
    trace_E = go.Scatter(
        name=f"{name} Decohesion Energy",
        x=df['separation distance (A)'],
        y=df['Decohesion energy (J/m^2)'],
        mode='lines+markers',
        yaxis='y1'
    )

    trace_S = go.Scatter(
        name=f"{name} Decohesion Stress",
        x=df['separation distance (A)'],
        y=df['Decohesion stress (GPa)'],
        mode='lines+markers',
        yaxis='y2'
    )
    layout = go.Layout(
        title=dict(
            text='Decohesion Energy and Stress',
            x=0.5,  # 标题居中
            xanchor='center'
        ),
        xaxis=dict(
            title_text="separation distance (A)",
            title_font=dict(
                size=18,
                color="#7f7f7f"
            )
        ),
        yaxis=dict(
            title="Decohesion energy (J/m^2)",
            titlefont=dict(
                size=18,
                color="#7f7f7f"
            )
        ),
        yaxis2=dict(
            title="Decohesion stress (GPa)",
            titlefont=dict(
                size=18,
                color="#7f7f7f"
            ),
            overlaying='y',
            side='right'
        )
    )
    trace = [trace_E, trace_S]
    fig = go.Figure(data=trace, layout=layout)
    fig.show()

    return [trace], layout


def dash_table(res_data: dict, decimal: int = 3, **kwargs) -> dash_table.DataTable:
    decohesion_e = [values[0] for values in res_data.values()]
    stress = [values[1] for values in res_data.values()]
    vacuum_size = [values[2] for values in res_data.values()]
    vacuum_size = [str(item) for item in vacuum_size]
    df = pd.DataFrame({
        "separation distance (A)": vacuum_size,
        "Decohesion energy (J/m^2)": round_format(decohesion_e, decimal),
        "Decohesion stress (GPa)": round_format([s / 1e9 for s in stress], decimal),
    })
    table = dash_table.DataTable(
        data=df.to_dict('records'),
        columns=[{'name': i, 'id': i} for i in df.columns],
        style_table={'width': TABLE_WIDTH,
                     'minWidth': TABLE_MIN_WIDTH,
                     'overflowX': 'auto'},
        style_cell={'textAlign': 'left'}
    )
    return table, df

