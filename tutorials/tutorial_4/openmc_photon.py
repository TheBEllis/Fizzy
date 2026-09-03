#!/usr/env/python3

import openmc
import openmc.lib
import os
import sys
import csv
import h5py


def make_model(compiled_source, mesh_id, temperature=296):

    model = openmc.Model()
    model.settings.dagmc = True
    # steel
    steel = openmc.Material(name="steel", temperature=temperature)

    steel.set_density("g/cm3", 8.0)
    steel.add_nuclide("B10", 0.000693, "wo")
    steel.add_nuclide("B11", 0.002807, "wo")
    steel.add_nuclide("C12", 0.039520273803427465, "wo")
    steel.add_nuclide("C13", 0.000479726196572539, "wo")
    steel.add_nuclide("Si28", 0.47, "wo")
    steel.add_nuclide("V50", 0.0003921925373160739, "wo")
    steel.add_nuclide("V51", 0.15960780746268394, "wo")
    steel.add_nuclide("Cr50", 0.73, "wo")
    steel.add_nuclide("Cr52", 14.07, "wo")
    steel.add_nuclide("Cr53", 1.6, "wo")
    steel.add_nuclide("Cr54", 0.4, "wo")
    steel.add_nuclide("Mn55", 1.14, "wo")
    steel.add_nuclide("Fe54", 3.95, "wo")
    steel.add_nuclide("Fe56", 62.51, "wo")
    steel.add_nuclide("Fe57", 1.46, "wo")
    steel.add_nuclide("Fe58", 0.2, "wo")
    steel.add_nuclide("Co59", 0.14, "wo")
    steel.add_nuclide("Ni58", 7.31, "wo")
    steel.add_nuclide("Ni60", 2.79, "wo")
    steel.add_nuclide("Ni61", 0.12, "wo")
    steel.add_nuclide("Ni62", 0.38, "wo")
    steel.add_nuclide("Ni64", 0.1, "wo")
    steel.add_nuclide("Mo92", 0.2974760972643644, "wo")
    steel.add_nuclide("Mo94", 0.1906159496693566, "wo")
    steel.add_nuclide("Mo95", 0.33284984504744625, "wo")
    steel.add_nuclide("Mo96", 0.3533051306488447, "wo")
    steel.add_nuclide("Mo97", 0.20516502475964524, "wo")
    steel.add_nuclide("Mo98", 0.5254922939545649, "wo")
    steel.add_nuclide("Mo100", 0.21509565865577795, "wo")
    steel.add_nuclide("Cu63", 0.06, "wo")
    steel.add_nuclide("Cu65", 0.03, "wo")

    materials = openmc.Materials()
    materials.append(steel)

    model.materials = materials

    surf = openmc.Sphere(r=250, surface_id=99999, boundary_type="vacuum")
    dag = openmc.DAGMCUniverse("/Projects/Fizzy/geometry/cube_dag.h5m", name='daguni', universe_id=42)

    cell = openmc.Cell(fill=dag, region=-surf, cell_id=99999)

    mesh = openmc.UnstructuredMesh(
        filename='/Projects/Fizzy/geometry/cube.e', library='libmesh', mesh_id=1)

    mesh_filter = openmc.MeshFilter(mesh)
    photon_filter = openmc.ParticleFilter("photon")

    photon_flux_tally = openmc.Tally(name="photon_flux")
    photon_flux_tally.filters = [mesh_filter, photon_filter]
    photon_flux_tally.scores = ['flux']
    photon_flux_tally.estimator = 'collision'

    model.tallies = [photon_flux_tally]

    model.geometry = openmc.Geometry([cell])
    model.settings.batches = 10
    model.settings.inactive = 0
    model.settings.output = {"summary": False,
                             "tallies": False}
    # 1e10 histories per batch
    model.settings.particles = 50000
    model.settings.photon_transport = True
    model.settings.electron_treatment = 'ttb'
    model.settings.run_mode = 'fixed source'
    model.settings.energy_mode = 'continuous-energy'
    # model.settings.sourcepoint = sourcepoint
    model.settings.max_tracks = 10000

    model.settings.source = openmc.CompiledSource(
        compiled_source, mesh_id)
    return model


model = make_model(
    '/Projects/FizzyCompiledSource/build/libCompiledSource.so', '1')

model.export_to_xml("./openmc_inputs.xml")
