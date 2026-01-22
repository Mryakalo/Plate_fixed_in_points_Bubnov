import numpy as np
import plotly.graph_objects as go

from app.helpers import ask_parameter


def plot_graph(w_function, a, b):
    # 1. Определяем количество точек разбиения вдоль x
    name = 'Number internal points between 0 and y'
    default = 11
    number_x_points = ask_parameter(name=name, params_type=int, default=default) + 2

    # 2. Определяем количество точек разбиения вдоль y
    name = 'Number internal points between 0 and y'
    default = 11
    number_y_points = ask_parameter(name=name, params_type=int, default=default) + 2

    # 3. Собирается массив со всеми перемещениями всех узлов
    matrix_w: np.ndarray = np.zeros((number_x_points, number_y_points))
    for i in range(number_x_points):
        for j in range(number_y_points):
            matrix_w[i][j] = w_function(a * i / (number_y_points - 1), b * j / (number_y_points - 1))
    # print(matrix_w)

    # 3. Заполняем массивы с координатами вдоль x и y
    x_axis = np.linspace(0, a, int(number_y_points))
    y_axis = np.linspace(0, b, int(number_y_points))
    z_axis = matrix_w * 1000

    y_axis, x_axis = np.meshgrid(y_axis, x_axis)

    # -------------------------------------
    # Create the 3D surface plot
    fig = go.Figure(data=[go.Surface(z=z_axis,
                                     x=x_axis,
                                     y=y_axis,
                                     colorscale='Viridis',
                                     colorbar_title_text=f"W, мм (max: {z_axis.max():.2f})"
                                     )])

    # Update layout
    fig.update_layout(
        title='Прогибы плиты',
        scene=dict(
            aspectmode='manual',
            aspectratio=dict(x=a/a, y=b/a, z=a/4/a),
            xaxis_title='x, м',
            yaxis_title='y, м',
            zaxis_title='W, мм',
            zaxis=dict(autorange="reversed",
                       visible=False),
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),

        ),
        font=dict(size=30, family="Arial", color='black'),
        autosize=False,
        width=1500,
        height=800,

    )

    fig.update_traces(contours={
        "z": {
            "show": True,         # Display z-axis contours
            "color": "white",     # Color of the contour lines
            "highlight": True,    # Highlight lines on hover
            "project": {"z": True}, # Project the lines onto the z-plane
        },
    }
        )


    fig.show()

