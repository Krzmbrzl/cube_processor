#!/usr/bin/env python3

from typing import List

import argparse
import numpy as np


class Atom:
    def __init__(self, atomic_number: int, nuclear_charge: float, position: np.ndarray) -> None:
        self.atomic_number = atomic_number
        self.nuclear_charge = nuclear_charge
        self.position = position


class CubeFile:
    def __init__(self, file_path: str) -> None:
        self.load_from(file_path)


    def load_from(self, file_path: str) -> None:
        """Parses the cube file whose path is given as the argument to this function"""
        with open(file_path, "r") as input_file:
            # Format spec from https://h5cube-spec.readthedocs.io/en/latest/cubeformat.html

            # First two lines are comments
            self.comment1 = input_file.readline()
            self.comment2 = input_file.readline()

            # {NATOMS (int)} {ORIGIN (3x float)} {NVAL (int)}?
            parts = input_file.readline().split()
            assert len(parts) in [4,5]

            nAtoms = abs(int(parts[0]))
            assert int(parts[0]) >= 0, "DSET_IDs not (yet) supported"

            self.origin = np.ndarray(shape=(3), dtype=float)
            for i in range(3):
                self.origin[i] = float(parts[i + 1])

            self.data_points_per_grid_point = int(parts[4]) if len(parts) > 4 else 1


            # {XAXIS (int) (3x float)}
            parts = input_file.readline().split()
            assert len(parts) == 4
            self.x_grid_points = int(parts[0])
            self.x_axis = np.ndarray(shape=(3), dtype=float)
            for i in range(3):
                self.x_axis[i] = float(parts[i + 1])
            # {YAXIS (int) (3x float)}
            parts = input_file.readline().split()
            assert len(parts) == 4
            self.y_grid_points = int(parts[0])
            self.y_axis = np.ndarray(shape=(3), dtype=float)
            for i in range(3):
                self.y_axis[i] = float(parts[i + 1])
            # {ZAXIS (int) (3x float)}
            parts = input_file.readline().split()
            assert len(parts) == 4
            self.z_grid_points = int(parts[0])
            self.z_axis = np.ndarray(shape=(3), dtype=float)
            for i in range(3):
                self.z_axis[i] = float(parts[i + 1])


            # Molecular geometry
            self.atoms: List[Atom] = []
            for i in range(nAtoms):
                # {GEOM (int) (float) (3x float)}
                parts = input_file.readline().split()
                assert len(parts) == 5

                self.atoms.append(Atom(atomic_number=int(parts[0]), nuclear_charge=float(parts[1]),
                                       position=np.array([float(parts[2]), float(parts[3]), float(parts[4])])))
            
            # {DSET_IDS (#x int)}
            # -> not yet supported
            
            # Read data
            # The data points are given as a sequence of whitespace-separated floating point numbers (written in scientific notation)
            # The values are given in the following hierarchy
            # - nPoints: All data points for a fixed grid point (x,y,z)
            # - zLine: All nPoints along the z-axis for fixed (x,y), starting at x=0, y=0
            # - yPlane: All zLines along the y-axis for fixed (x) starting at x=0
            # - xCube: All yPlanes along the x-axis
            self.data: np.ndarray = np.fromstring(input_file.read(), dtype=float, sep=' ')

            assert len(self.data) == self.x_grid_points * self.y_grid_points * self.z_grid_points * self.data_points_per_grid_point


    def min_value(self) -> float:
        """Gets the minimum value in this data set"""
        return np.min(self.data)

    def abs_min_value(self) -> float:
        """Gets the data value in this data set that is closest to zero"""
        return np.min(np.abs(self.data))

    def max_value(self) -> float:
        """Gets the maximum value in this data set"""
        return np.max(self.data)

    def abs_max_value(self) -> float:
        """Gets the data value in this data set that is furthest away zero"""
        return np.max(np.abs(self.data))

    def summed_data(self) -> float:
        """Gets the sum over all data points"""
        return np.sum(self.data)

    def abs_summed_data(self) -> float:
        """Gets the sum over all absolute values of the data points"""
        return np.sum(self.data)

    def volume(self) -> float:
        """Gets the volume of the cuboid represented by this object. The unit of the returned volume is the third power of the length unit used in the
        underlying cube file, which (according to the specs) should be Bohr."""
        return np.linalg.norm(self.x_axis) * np.linalg.norm(self.y_axis) * np.linalg.norm(self.z_axis)

    def integrate(self) -> float:
        """Gets the property represented by this cube file, integrated over the entire space of the cuboid represented here. The unit of the returned
        quantity depends on the unit of whatever data is represented in this cube file. It is that times the cube of the length units used in this
        file (normally Bohr)"""
        return self.summed_data() * self.volume()

    def isosurface_threshold_value(self, coverage_percent: float = 90) -> float:
        """Gets the threshold value at which an isosurface has to be drawn, such that the volume enclosed inside said isosurface (potentially bounded
        by the boundaries of this cuboid) encloses a given percentage of the total property represented by this cube file"""
        if coverage_percent >= 100:
            return 0

        sorted_indices = np.argsort(self.data)

        target_sum = self.summed_data() * coverage_percent / 100

        current_sum = 0
        for i in reversed(sorted_indices):
            current_sum += self.data[i]

            if current_sum >= target_sum:
                # The threshold must be chosen from the original, unmodified data set
                return self.data[i]

        return self.data[sorted_indices[-1]]


def main() -> None:
    parser = argparse.ArgumentParser(description="CLI tool for processing Gaussian cube files")
    
    parser.add_argument("cube_file", help="The cube file to process")
    parser.add_argument("--print-info", help="Print some general info about the read cube file", action="store_true", default=False)
    parser.add_argument("--calc-iso-value",
        help="Determine a threshold value for a n isosurface such that the enclosed volume will contain the given percentage of the total property",
        const=90, nargs="?", type=int, metavar="PERCENTAGE")

    args = parser.parse_args()

    cube = CubeFile(file_path = args.cube_file)

    if args.print_info:
        print("Summary for cube file '%s':" % args.cube_file)
        print(" - Min. value:      {:+.2e}".format(cube.min_value()))
        print(" - Abs. min. value: {:+.2e}".format(cube.abs_min_value()))
        print(" - Max. value:      {:+.2e}".format(cube.max_value()))
        print(" - Abs. max. value: {:+.2e}".format(cube.abs_max_value()))
        print(" - Sum(data):       {:+.4e}".format(cube.summed_data()))
        print(" - Integrated data: {:+.4e}".format(cube.integrate()))

    if args.calc_iso_value:
        print("To create an isosurface enclosing {:d}% of the contained property, use a threshold of\n{:.4e}".format(
                args.calc_iso_value,
                cube.isosurface_threshold_value(coverage_percent=args.calc_iso_value)
            )
        )


if __name__ == "__main__":
    main()
