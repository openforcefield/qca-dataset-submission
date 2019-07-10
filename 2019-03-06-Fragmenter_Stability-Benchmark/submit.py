import qcfractal.interface as ptl
import json
import tarfile

# Build a client to the target server
client = ptl.FractalClient.from_file("localhost:7777")

# Unpack the data
with tarfile.open('stability_benchmark_inputs.tar.gz', 'r') as f:
    td_inputs = json.load(f.extractfile('stability_benchmark_inputs.json'))

# Format from raw data to a TorsionDriveDataset
ds = ptl.collections.TorsionDriveDataset("Fragment Stability Benchmark", client=client)
for td in td_inputs:
    label = td
    qcmol = td_inputs[td]['initial_molecule']
    dih = td_inputs[td]['dihedral']
    grid = td_inputs[td]['grid']
    attributes = td_inputs[td]['identifiers']
    ds.add_entry(name=label, initial_molecules=qcmol, dihedrals=dih, grid_spacing=grid, attributes=attributes)

# Build out specification of compute
optimization_spec = {
    'program': 'geometric',
    'keywords': {
       'coordsys': 'tric',
       'enforce': 0.1,
       'reset': True,
       'qccnv': True,
       'epsilon': 0
         }
      }

qc_spec = {'driver': "gradient",
  'method': 'b3lyp-d3(bj)',
  'program': 'psi4',
  'basis': 'dzvp',
  'keywords': 2} # Keywords 2 id map compute wiberg bond orders, dipoles, and quadrupoles.

# Examine keywords with
# client.query_keywords(id=2)[0].dict()

# Add the specification under the "default" name
ds.add_specification("default", optimization_spec, qc_spec, description="Default OpenFF B3LYP-D3(BJ)/DZVP values", overwrite=True)

# Add the computations
ds.compute("default", tag="openff")
