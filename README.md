# Decomposition of large-scale CVRP Instances

This repository presents a proof-of-concept for including instance decomposition within two popular metaheuristics for the Capacitated Vehicle Routing Problem (CVRP):
* The Adaptive Large Neighbourhood Search of Stefan Ropke
* The Hybrid Genetic Algorithm of Thibaut Vidal

We propose several decomposition techniques, which we describe in detail in our preprint:

```
@techreport{cvrp_decomposition,
    title={Decomposition strategies for vehicle routing heuristics},
    author={Santini, Alberto and Schneider, Michael and Vidal, Thibaut and Vigo, Daniele},
    year=2022,
    url={https://santini.in/files/papers/santini-schneider-vidal-vigo-2022.pdf}
}
```

## Authors

The original ALNS code is by Stefan Ropke.
Stefan wrote most of the code contained in directory `src/alns`.
Alberto Santini then edited Stefan's code and implemented instance decomposition.
The original paper in which Stefan introduced ALNS is the following:

```
@article{Ropke2006a,
    title={A unified heuristic for a large class of vehicle routing problems with backhauls},
    author={Ropke, Stefan and Pisinger, David},
    year=2006,
    journal={European Journal of Operational Research},
    volume=171,
    number=3,
    pages={750--775},
    doi={10.1016/j.ejor.2004.09.004}
}
```

The original HGS code is by Thibaut Vidal.
Thibaut also released his code in repository [vidalt/HGS-CVRP](https://github.com/vidalt/HGS-CVRP) and through the below paper.
Alberto Santini then edited Thibaut's code and implemented instance decomposition.

```
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

```
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