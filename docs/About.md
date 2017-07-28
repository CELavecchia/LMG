---
layout: page
title: "About"
logo: "img/home-bg.jpg"
description: "About LMG"
header-img: "img/home-bg.jpg"
---

![LMG_logo](img/logos/logo.png){:width="100%"}
Cite LMG: <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.835756.svg" alt="DOI_logo" width="50%">

# Table of contents
- [Project Summary](#Summary)  
- [Installation](#Installation)  
- [Getting started](#Start)  
- [License](#License)  
- [Contributing](#Contributing)  
- [Code of conduct](#CodeOfConduct)  
- [Road Map](#RoadMap)  

## Project summary <a name="Summary"></a>
Lumbar Model Generator (LMG) is a MATLAB toolbox for semi-automatic generation of lumbar finite element geometries.
This toolbox allows to obtain:
- the geometrical model of the lumbar spine (from the vertebrae L1 to the L5 including the intervertebral disc IVD);
- the surface models of the bodies involved (STL files);
- the solid meshed model, generated with hexahedral elements for the IVD and tetrahedral elements for the vertebrae;

A work in progress feature is the pre-processing of the solid meshed model to prepare the .feb file to run the finite element simulation in FEBio.

![LMG overview](img/overview.png){:width="100%"}

# Installation <a name="Installation"></a>  
### 1. Installing 3rd party packages
Skip this step if finite element analysis with FEBio is not required.

| Package | Purpose | Included? | Download |
|:--|:--|:--:|--:|
|[__FEBio__](https://febio.org) <br/> [![FEBio](/img/logos/febioLogo.png){:height="100px"}](https://febio.org)|FEBio is a finite element solver and is used in LMG and GIBBON for all finite element analysis. |__No__|[__FEBio website__](https://febio.org) |
|<br/> [__GIBBON__](https://gibboncode.org) <br/> [![export_fig](/img/logos/gibbonlogo.png){:height="100px"}](https://gibboncode.org)| <br/> GIBBON The Geometry and Image-Based Bioengineering add-On |__No__|[__GIBBON website__](https://gibboncode.org) |
|<br/> [__TetGen__]() <br/> [![tetGen](/img/logos/tetgenLogo.gif){:height="100px"}](https://wias-berlin.de/software/tetgen/)| <br/> Is used for tetrahedral meshing (and possibly constrained 3D Delaunay tessellation). See for instance `HELP_runTetGen.m`|__Yes__| For other versions: [__TetGen website__](https://wias-berlin.de/software/tetgen/)|

### 2. Run `install_LMG.m`
By running `install_LMG.m` the LMG and FEBio (if needed) path definitions will be added and saved to MATLAB. The help and documentation will also be integrated. Once finished you will be asked to __restart MATLAB__. `install_LMG.m` can be found in the main LMG folder.

## Getting started <a name="Start"></a>

##### Access the integrated help
* To access the help documentation from MATLAB click on the HELP browser then click on `LMG toolbox` under `Supplemental Software`. This will open the toolbox help and documentation which is now searchable and integrated just like the rest of MATLAB's help and documentation.  

## License <a name="License"></a>
[![License](https://img.shields.io/badge/License-BSD%203--Clause-green.svg)](https://github.com/CELavecchia/LMG/blob/master/LICENSE)

## Contributing <a name="Contributing"></a>
Coming soon

## Code of conduct <a name="CodeOfConduct"></a>
Coming soon

## Roadmap <a name="RoadMap"></a>
Coming soon
