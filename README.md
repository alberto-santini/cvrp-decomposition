# Decomposition Strategies for Vehicle Routing Heuristics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6097671.svg)](https://doi.org/10.5281/zenodo.6097671)

This repository presents a proof-of-concept for including instance decomposition within two popular metaheuristics for the Capacitated Vehicle Routing Problem (CVRP):
* The Adaptive Large Neighbourhood Search of Stefan Ropke
* The Hybrid Genetic Algorithm of Thibaut Vidal

We propose several decomposition techniques, which we describe in detail in our paper:

```bib
@article{cvrp_decomposition,
    title={Decomposition strategies for vehicle routing heuristics},
    author={Santini, Alberto and Schneider, Michael and Vidal, Thibaut and Vigo, Daniele},
    year=2023,
    doi={10.1287/ijoc.2023.1288},
    journal={{INFORMS Journal on Computing}},
    volume=35,
    number=3,
    pages={543--559}
}
```

You can cite this repository itself through Zenodo:

```bib
@article{cvrp_decomposition_2022,
    title={Repository cvrp-decomposition},
    author={Alberto Santini},
    year=2022,
    doi={10.5281/zenodo.6097671},
    url={https://github.com/alberto-santini/cvrp-decomposition}
}
```

## Authors

The original ALNS code is by Stefan Ropke.
Stefan wrote most of the code contained in directory `src/alns`.
Alberto Santini then edited Stefan's code and implemented instance decomposition.
The original paper in which Stefan introduced ALNS is the following:

```bib
@article{ropke_2006,
    title={An adaptive large neighborhood search heuristic for the pickup and delivery problem with time windows},
    author={Ropke, Stefan and Pisinger, David},
    year=2006,
    journal={Transportation Science},
    volume=40,
    number=4,
    pages={455--472},
    doi={10.1287/trsc.1050.0135}
}
```

The original HGS code is by Thibaut Vidal.
Thibaut also released his code in repository [vidalt/HGS-CVRP](https://github.com/vidalt/HGS-CVRP) and through the below paper.
Alberto Santini then edited Thibaut's code and implemented instance decomposition.

```bib
@article{vidal_2022,
    title={Hybrid genetic search for the CVRP: Open-source implementation and SWAP\textsupersctipt{*} neighborhood},
    author={Vidal, Thibaut},
    year=2022,
    journal={Computers \& Operations Research},
    volume=140,
    doi={10.1016/j.cor.2021.105643}
}
```

## Usage

The preferred way to build the two programmes (one for ALNS and one for HGS) is trhough cmake.
We provide a `CMakeList.txt` file, which allows to build the executables as follows:

```
$ git clone https://github.com:alberto-santini/cvrp-decomposition.git
$ cd cvrp-decomposition
$ mkdir build
$ cd build
$ cmake -DCMAKE_BUILD_TYPE=Release ..
$ cmake
$ make
```

## Instances

Instances contained in folder `data` are part of the so-called "Uchoa dataset".
We refer to the following paper for further information on their generation and characteristics:

```bib
@article{uchoa_2017,
    title={New benchmark instances for the capacitated vehicle routing problem},
    author={Uchoa, Eduardo and Pecin, Diego and Pessoa, Artur and Poggi, Marcus and Vidal, Thibaut and Subramanian, Anand},
    year=2017,
    journal={European Journal of Operational Research},
    volume=257,
    number=3,
    pages={845--858},
    doi={10.1016/j.ejor.2016.08.012}
}
```

## Acknowledgements

We are extremely grateful to Stefan Ropke for sharing with us his code for ALNS applied to the CVRP.

## License

The original HGS-CVRP is released under the MIT license, which we report in file `LICENSE`.
