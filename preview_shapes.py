from argparse import ArgumentParser
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy.core.defchararray import startswith
from numpy.core.function_base import linspace
import pandas as pd
import numpy as np
from typing import Protocol
from dataclasses import dataclass, field
from enum import Enum, auto


class AxiCodendType(Enum):
    T0 = auto()
    T90 = auto()


@dataclass
class AxiCodendResult():
    meshes_along: int
    meshes_across: int
    codend_type: 0
    meridian: list
    n_nodes: int = field(init=False)
    n_dof: int = field(init=False)
    
    def __post_init__(self):
        self.n_nodes = 1 + 4 * self.meshes_along
        self.n_dof = self.n_nodes * 2 
        if self.codend_type == AxiCodendType.T90:
            self.n_nodes = 1 * 2* self.meshes_along
            self.n_dof = self.n_nodes * 3


class CodendDrawer(Protocol):
    def dof_to_coords(self, obj: AxiCodendResult) -> np.ndarray:
        ...

    def draw(self, obj: AxiCodendResult):
        ...


class T0Drawer:
    def dof_to_coords(self, codend: AxiCodendResult, theta: float = 0, clockwise: bool = True) -> np.ndarray:
        """
        Expand data and transform meridian from axis-symetric coordinates to 3d
        [DOF, 1, 2] * [DOF, 2, 3] = [DOF, 1, 3]
        """

        n_nodes = 1 + 4 * codend.meshes_along 
        T_start = np.array([
            [1, 0, 0],
            [0, np.sin(theta), np.cos(theta)]
        ])

        theta_inc = np.pi / codend.meshes_across
        T_end = np.array([
            [1, 0, 0],
            [0, np.sin(theta + theta_inc), np.cos(theta + theta_inc)]
        ])
        if not clockwise:
            T_start, T_end = T_end, T_start

        T = np.zeros((n_nodes, 2, 3))
        T[0,:,:] = T_start
        T[1:] = np.tile(
            np.stack([T_start, T_end, T_end, T_start], axis=0),
            (codend.meshes_along, 1, 1))

        data2d = np.reshape(codend.meridian , (-1, 2))
        batch_data = np.expand_dims(data2d, axis=1)
        return np.matmul(batch_data, T).squeeze(); 


    def draw(self, codend: AxiCodendResult, draw_half=True):
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        n_meridians = int(0.5 * codend.meshes_across) if draw_half else codend.meshes_across
        data3d_all = np.zeros(((codend.n_nodes + 1) * n_meridians, 3))
        for i, theta in enumerate(np.linspace(0, np.pi, n_meridians)):
            line_data = self.dof_to_coords(codend, theta)
            ax.plot(line_data[:,0], line_data[:,1], line_data[:,2], 'g', linewidth=1)
    
        set_axes_equal(ax)
        plt.show()


def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])


def parse_shapes_txt(path: str, sep="\t") -> pd.DataFrame:
    df = pd.read_csv(path, sep, header=None)

    results = []
    for column in df:
        meridian = df[column].to_numpy()
        results.append(AxiCodendResult(100, 100, AxiCodendType.T0, meridian))
    return results


def main(shapes_file: str, results_file: str=None):
    results = parse_shapes_txt(shapes_file)
    drawer = T0Drawer()
    drawer.draw(results[1])

    # dc = DrawController(drawer, results)
    # dc.draw() will get an interface with arrows

if __name__ == "__main__":
    ap = ArgumentParser("Preview for AxiCodend output shapes and results")
    ap.add_argument("-s", "--shapes-file", required=True, type=str,
        help="Path to the file containing resulting shapes.")
    ap.add_argument("-r", "--results-file", required=False, type=str,
        help="Path to the file containing summary results.")
    args = ap.parse_args()

    main(args.shapes_file, args.results_file)