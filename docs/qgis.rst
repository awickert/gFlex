QGIS Processing Provider
========================

`processing_gflex <https://github.com/awickert/processing_gflex>`_ exposes
gFlex as a no-code algorithm in the QGIS Processing Toolbox and Graphical
Modeler.  If your load and elastic-thickness data are already rasters in a
QGIS project, this is the easiest way to run a flexural calculation without
writing any Python.

**Requirements**

* QGIS ≥ 3.16
* gFlex ≥ 2.0.0 (installed automatically on first use if not already present)

**Capabilities**

* All 2-D solution methods: ``fd``, ``fft``, ``sas``
* Scalar or raster (spatially variable) :math:`T_e`
* All boundary conditions (``infinite``, ``clamped``, ``pinned``, ``free``,
  ``mirror``, ``periodic``)
* In-plane stresses (:math:`\sigma_{xx}`, :math:`\sigma_{yy}`,
  :math:`\sigma_{xy}`)
* Raster output: lithospheric deflection [m]

**Installation**

Install processing_gflex from the QGIS Plugin Manager (search
"processing_gflex") or from the command line::

   pip install "gflex>=2.0.0"   # prerequisite; also installed automatically on first use

After enabling the plugin, the *gFlex Flexure* algorithm appears in the
Processing Toolbox under **gFlex**.

**Usage**

Open the Processing Toolbox (**Processing → Toolbox**), locate
**gFlex → Compute lithospheric flexure**, and fill in the dialog:

* **Load raster** — surface-normal stress [Pa], :math:`\rho_\text{load}\,g\,h`
* **Elastic thickness** — scalar value [m] or raster (variable :math:`T_e`)
* **Method** — ``fd``, ``fft``, or ``sas``
* **Boundary conditions** — one per edge (north, south, east, west)
* **Physical parameters** — :math:`E`, :math:`\nu`, :math:`g`,
  :math:`\rho_m`, :math:`\rho_\text{fill}`
* **Output deflection** — path for the output raster

The same algorithm is accessible headlessly via ``processing.run()``
in the QGIS Python console or in a Graphical Modeler chain::

   import processing
   result = processing.run(
       "gflex:computeflexure",
       {
           "INPUT_LOAD": "/path/to/load.tif",
           "INPUT_TE": 35000,           # scalar 35 km
           "METHOD": "fd",
           "BC_WEST": "infinite",
           "BC_EAST": "infinite",
           "BC_NORTH": "infinite",
           "BC_SOUTH": "infinite",
           "OUTPUT": "/path/to/deflection.tif",
       },
   )

See the `processing_gflex repository <https://github.com/awickert/processing_gflex>`_
for full parameter documentation and changelog.
