# SeDuS - Segmental Duplication Simulator

## Overview

SeDuS (Segmental Duplication Simulator) is the first flexible and user-friendly forward-in-time simulator designed to model patterns of molecular evolution within segmental duplications that undergo interlocus gene conversion and crossover. SeDuS introduces known features of interlocus gene conversion, including biased directionality and dependence on local sequence identity. Additionally, the simulator incorporates various aspects such as different selective pressures acting upon copy number and flexible crossover distributions.

The simulator features a graphical user interface (GUI) that allows users to fine-tune relevant parameters quickly and perform real-time analysis of the evolution of duplicated segments.

## Features

- **Interlocus Gene Conversion Simulation**: Models biased directionality and local sequence identity dependence.
- **Crossover Simulation**: Incorporates flexible crossover distributions within segmental duplications.
- **Selective Pressure Modelling**: Simulates different selective pressures acting on copy number variation.
- **Graphical User Interface (GUI)**: Provides an intuitive interface for parameter adjustment and real-time analysis.
- **Cross-Platform Availability**: Available for Linux, OS X, and Windows.

## Availability and Implementation

SeDuS is implemented in C++ and offers both command line and graphical user interface (GUI) options. The GUI is developed using Qt C++ to ensure a user-friendly experience across multiple platforms.

- **Source Code and Binaries**: Source code and precompiled binaries for Linux, OS X, and Windows are freely available at [www.biologiaevolutiva.org/sedus/](http://www.biologiaevolutiva.org/sedus/).
- **Tutorial**: A comprehensive tutorial with detailed descriptions of the implementation, parameters, and output files is available online.

## Installation

### Prerequisites

- **C++ Compiler**: Required if you choose to compile the source code yourself.
- **Qt Libraries**: Necessary for running the graphical user interface (only if compiling from source).
  
### Building from Source

1. **Clone** the repository:
   ```bash
   git clone https://github.com/jmanuelfuentes/sedus.git
   cd sedus
