# MSAT
MSAT aims to provide a helpful tool to quantitatively analyze the shock instability problem for the second-order finite-volume MUSCL schemes. Detailed information about the matrix stability analysis method for the first-order scheme can be found in the following: 
 - [M. Dumbser, J.-M. Moschetta, J. Gressier, A matrix stability analysis of the carbuncle phenomenon, Journal of Computational Physics 197 (2) (2004) 647–670](https://doi.org/10.1016/j.jcp.2003.12.013)

And in this software, we extend the matrix stability analysis method for the second-order MUSCL scheme.

## Authors
- Weijie Ren <rwj@nudt.edu.cn>
- Wenjia Xie <xiewenjia@nudt.edu.cn>
- Ye Zhang <zhangye17a@nudt.edu.cn>
- Hang Yu <yuhang18@gfkd.edu.cn>
- Zhengyu Tian <tianzhengyu_kd@163.com>

## Usage
Following the next step to use MSAT:
 1. Clone the directory with <https://github.com/JameRwj/MSAT.git>
 2. Change the settings in *Settings.dat*
 3. Generate the computational grid and export to *Grid.dat*
 4. Execute *Main Program* to run the software

Note that MSAT is written in Fortran and uses MKL to calculate the eigenvalues and eigenvectors of the stability matrix, so the Fortran compiler and API should be prepared in advance. The software is tested on Windows with Visual Studio 2017, Intel Fortran compiler 2021, and Intel oneAPI 2022.
