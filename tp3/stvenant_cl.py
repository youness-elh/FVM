#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""resolution of a transport equation by the finite volume method on regular grid"""


from __future__ import absolute_import, print_function
import pyopencl as cl
import numpy as np

# Load common utilities
import sys
sys.path.append('..')
from utils import Figure, load_kernel, get_ite_title, parse_cl_args


## Default values

# number of conservative variables
_m = 3

# grid size
_nxy = 1024
_nx = _nxy
_ny = _nxy
_Lx = 1.
_Ly = 1.

gpes = 9.81
umax = 1.5
hmax = 2.

vmax = np.sqrt(gpes * hmax) + umax

# time stepping
cfl = 0.8
_Tmax = 0.05


def solve_ocl(m=_m, nx=_nx, ny=_ny, Tmax=_Tmax, Lx=_Lx, Ly=_Ly, animate=True, ivplot=0, precision="single", **kwargs):

    dx = Lx / nx
    dy = Ly / ny
    dt = cfl * (dx * dy) / (2 * dx + 2 * dy) / vmax

    # For plotting
    x = np.linspace(0., Lx, num=nx)
    y = np.linspace(0., Ly, num=ny)

    # load and adjust  C program
    parameters = {'nx': nx,
                  'ny': ny,
                  'dx': dx,
                  'dy': dy,
                  'dt': dt,
                  'm': m,
                  'lambda': vmax/2.,
                  }

    np_real, source = load_kernel("stvenant_kernels.cl", parameters, precision=precision, print_source=True,
                                  module_file=__file__)

    # OpenCL init
    ctx = cl.create_some_context()
    mf = cl.mem_flags

    # compile OpenCL C program
    prg = cl.Program(ctx, source).build(options="-cl-fast-relaxed-math")

    # create OpenCL buffers
    buffer_size = m * nx * ny * np.dtype(np_real).itemsize
    wn_gpu = cl.Buffer(ctx, mf.READ_WRITE, size=buffer_size)
    wnp1_gpu = cl.Buffer(ctx, mf.READ_WRITE, size=buffer_size)

    # create a queue (for submitting opencl operations)
    queue = cl.CommandQueue(ctx, properties=cl.command_queue_properties.PROFILING_ENABLE)

    # init data
    event = prg.init_sol(queue, (nx * ny, ), (32, ), wn_gpu)
    event.wait() # on attend la fin de l'exec du kernel

    # number of animation frames
    nbplots = 120
    itemax = int(np.floor(Tmax / dt))
    iteplot = int(itemax / nbplots)

    #iteplot = 1
    # time loop
    t = 0
    ite = 0
    elapsed = 0.
    wn_cpu = np.empty((m * nx * ny, ), dtype=np_real)

    if animate:
        fig = Figure(title=r"$n_x = {}, n_y = {}$".format(nx, ny),
                     levels=np.linspace(0.5, 2., 64))

    print("start OpenCL computations...")
    while t < Tmax:
        t += dt
        # queue : file d'attente pour réaliser les calculs
        # nx*ny = nb_proc opencl (1D)
        # 64 : batch size (nx*ny proc répartis en groupes de 64 proc)
        event = prg.time_step(queue, (nx * ny, ), (64, ), wn_gpu, wnp1_gpu)
        event.wait()
        elapsed += 1e-9 * (event.profile.end - event.profile.start)
        # exchange buffer references for avoiding a copy
        wn_gpu, wnp1_gpu = wnp1_gpu, wn_gpu
        ite_title = get_ite_title(ite, t, elapsed)

        if ite % iteplot == 0 and animate:
            cl.enqueue_copy(queue, wn_cpu, wn_gpu).wait()
            wplot = np.reshape(wn_cpu, (m, nx, ny))
            fig.update(x, y, wplot[ivplot, :, :], suptitle=ite_title, cb=ite == 0)
        else:
            print(ite_title, end='\r')
        ite += 1

    # copy OpenCL data to CPU and return the results
    cl.enqueue_copy(queue, wn_cpu, wn_gpu).wait()

    wplot_gpu = np.reshape(wn_cpu, (m, nx, ny))
    return x, y, wplot_gpu


if __name__ == '__main__':
    args = parse_cl_args(n=_nxy, L=1.0, tmax=_Tmax, description='Solve Saint-Venant equations using PyOpenCL')

    # gpu solve
    x, y, wplot_gpu = solve_ocl(**vars(args))

    fig = Figure(title="Final solution")
    ivplot = 0
    fig.update(x, y, wplot_gpu[ivplot, :, :], cb=True, show=True)
