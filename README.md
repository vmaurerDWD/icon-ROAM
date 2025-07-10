## icon-ROAM

This repository contains parts of the source code of the release [icon_2024.07](https://gitlab.dkrz.de/icon/icon-model/-/tree/release-2024.07-public), 
with modifications as in https://gitlab.dkrz.de/icon/icon-nwp/-/tree/icon-ROAM (access to the icon-nwp repository is only possible with DKRZ credentials and access rights to https://gitlab.dkrz.de/icon/icon-nwp).
The parts of the source code as made available in this github repository is containing interfaces to the OASIS3-MCT coupler based on 
[Hagemann et al. (2024)](https://gmd.copernicus.org/articles/17/7815/2024/gmd-17-7815-2024.html), which can be used to couple the ocean model NEMO.

The configure scripts are modified compared to the version of the ICON public release so that the coupling via OASIS is switched on by default (`--enable_coupling_OAS=yes`). 
To compile the code, download [icon_2024.07](https://gitlab.dkrz.de/icon/icon-model/-/tree/release-2024.07-public) and replace
- `src/`
- `config/`
- `configure.ac`

with the versions in this repository. For compilation with the OASIS,
a compiled version of OASIS3-MCTv5 and an adapted configure wrapper are necessary. See the `*_oasis` versions in `config/dwd` and `config/ecmwf` as working examples on how to add the OASIS library paths.
For common build instructions, see
https://gitlab.dkrz.de/icon/icon-model/-/blob/release-2024.07-public/doc/Quick_Start.md .


## Partnership

The ICON partnership coordinates research activities developing, maintaining, and supporting the ICON modeling framework. 
ICON partner institutions are:
- [DWD](https://www.dwd.de/EN/Home/home_node.html)
- [MPI-M](https://www.mpimet.mpg.de/en/home/)
- [DKRZ](https://www.dkrz.de/en/dkrz-partner-for-climate-research?set_language=en)
- [KIT](https://www.kit.edu/english/index.php)
- [C2SM](https://c2sm.ethz.ch/): [ETH](https://ethz.ch/en.html) and [MeteoSwiss](https://www.meteoswiss.admin.ch/)
More information about ICON is available in the [project's public web page](http://icon-model.org).

## License

ICON is available under a BSD 3-clause license. See [LICENSES/](./LICENSES) for license information and [AUTHORS.TXT](./AUTHORS.TXT) for a list of authors.
