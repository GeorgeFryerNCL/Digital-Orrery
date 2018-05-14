"""
Sets up an N body problem using initial conditions from a file, and
produces a graph and output file of the motion of the bodies

George Fryer 2018
"""


import matplotlib.animation as anim
import matplotlib.pyplot as pplot
import numpy
import scipy.constants
import scipy.integrate



class NBody:
    """
    Stores the bodies, and provides a
    method to compute the motion
    """


    days_per_year = 365
    # Use units of km^3.kg^-1.day^-2
    G = scipy.constants.G * (60*60*24)**2 * 10**(-9)


    def __init__(self, bodies):
        """
        Accepts as arguments a list of objects of the Body class, describing
        the initial conditions of the N bodies
        """

        self._bodies = bodies
        self._N = len(bodies)


    def run(self, time, dt, datafile="out", graphfile="out", videofile="out", force_video=False,
            ffmpeg_exec="C:\\Program Files\\FFmpeg\\bin\\ffmpeg.exe"):
        """
        Computes the motion of the bodies for an integer number of days specified
        by the time parameter, in steps of dt days. The graph
        and the data it represents are outputted to graphfile and datafile
        """

        N = self._N
        m = numpy.array([body.m for body in self._bodies])
        r = numpy.array([body.r for body in self._bodies])
        v = numpy.array([body.v for body in self._bodies])
        M = numpy.sum(m)

        # Move centre of mass to origin
        R = numpy.array([numpy.sum([m[i]*r[i][0] for i in range(N)])/M,
                         numpy.sum([m[i]*r[i][1] for i in range(N)])/M])
        for r_elem in r:
            r_elem[0] = r_elem[0] - R[0]
            r_elem[1] = r_elem[1] - R[1]

        # Calculate velocity of centre of mass
        v_R = numpy.array([numpy.sum([m[i]*v[i][0] for i in range(N)])/M,
                           numpy.sum([m[i]*v[i][1] for i in range(N)])/M])

        # Set up the ODEs. U = [r_1_x, v_1_x, r_1_y, v_1_y,
        #                       r_2_x, v_2_x, r_2_y, v_2_y,
        #                       r_3_x, v_3_x, r_3_y, v_3_y,
        #                       ...]
        def dU_dt(U, t):
            derivs = []
            for i in range(N):
                force_sum_x = 0; force_sum_y = 0
                for j in range(N):
                    if not i == j:
                        force_sum_x = force_sum_x + m[j]*(U[4*i]-U[4*j]) \
                                        / (((U[4*i]-U[4*j])**2 + (U[4*i + 2]-U[4*j + 2])**2)**(3/2))
                        force_sum_y = force_sum_y + m[j]*(U[4*i + 2]-U[4*j + 2]) \
                                        / (((U[4*i]-U[4*j])**2 + (U[4*i + 2]-U[4*j + 2])**2)**(3/2))
                derivs.append(U[4*i + 1] - v_R[0])
                derivs.append(-NBody.G*force_sum_x)
                derivs.append(U[4*i + 3] - v_R[1])
                derivs.append(-NBody.G*force_sum_y)
            return derivs

        # Solve the ODEs
        time_steps = [i*dt for i in range(time//dt + 1)]
        initial = []
        for i in range(N):
            initial.append(r[i][0])
            initial.append(v[i][0])
            initial.append(r[i][1])
            initial.append(v[i][1])
        soln = scipy.integrate.odeint(dU_dt, initial, time_steps, rtol=1e-10)

        # Use soln to calculate r, v, and KE
        r_soln = numpy.array([numpy.array(list(zip(soln[:,4*i],soln[:,4*i + 2]))) for i in range(N)])
        v_soln = numpy.array([numpy.array(list(zip(soln[:,4*i + 1],soln[:,4*i + 3]))) for i in range(N)])
        KE = []
        for i in range(len(time_steps)):
            KE.append(sum([(m[j]/2)*(v_soln[j][i][0]**2 + v_soln[j][i][1]**2) for j in range(N)]))

        # Output results to file
        self._write_txt(datafile, time_steps, r_soln, v_soln, KE)

        # Produce static graph
        self._write_graph(graphfile, time_steps, r_soln, KE)

        # Produce video animation only for short simulations <1k time steps
        vid_frame_limit = 1e3
        if force_video or len(time_steps) <= vid_frame_limit:
            self._write_video(videofile, ffmpeg_exec, time_steps, r_soln)


    def _write_txt(self, datafile, time_steps, r, v, KE):
        """
        Outputs data for all time steps to datafile
        """

        # Main data file
        with open(datafile + ".dat", mode="w") as out:
            if not out:
                raise IOError("Could not open output file")

            for i in range(self._N):
                out.write(self._bodies[i].name \
                    + "\nTIME (days)\t\tX (km)\t\tY (km)\t\tv_X (km/day)\t\tv_Y (km/day)\n")
                for j in range(len(time_steps)):
                    out.write(str(time_steps[j]) + "\t\t" +
                          str(r[i][j][0]) + "\t\t" +
                          str(r[i][j][1]) + "\t\t" +
                          str(v[i][j][0]) + "\t\t" +
                          str(v[i][j][1]) + "\n")

        # Kinetic energy data file
        with open(datafile + "_ke.dat", mode="w") as out:
            if not out:
                raise IOError("Could not open output file")
            out.write("\nTIME (days)\t\tKINETIC ENERGY (J)\n")
            for j in range(len(time_steps)):
                out.write(str(time_steps[j]) + "\t\t" + str(KE[j]) + "\n")


    def _write_graph(self, graphfile, time_steps, r, KE):
        """
        Outputs data for all time steps as a graph to graphfile and
        outputs graph of total kinetic energy
        """

        dpi_g = 1000
        dpi_ke_g = 200

        # Position graph
        pplot.axis("equal")
        for i in range(self._N):
            pplot.plot(r[i][:,0], r[i][:,1],
                       self._bodies[i].marker,
                       linewidth=0.4,
                       color=self._bodies[i].colour,
                       label=self._bodies[i].name)
        pplot.xlabel("X (km)"); pplot.ylabel("Y (km)")
        sim_len = time_steps[-1]
        time_str = str(sim_len) + " days" if sim_len < 10*NBody.days_per_year \
                    else ("~" if sim_len % NBody.days_per_year else "") \
                        + str(int(sim_len//NBody.days_per_year)) + " years"
        pplot.title(str(self._N) + "-body system over " + time_str)
        pplot.legend(loc=1)
        pplot.rcParams["agg.path.chunksize"] = 1000
        pplot.savefig(graphfile + ".png", dpi=dpi_g)

        # Total kinetic energy graph
        margin = 1.1
        pplot.clf()
        pplot.gca().set_ylim(0, margin*max(KE))
        time_steps,time_unit = (time_steps,"(days)") if sim_len < 10*NBody.days_per_year \
                                else ([t/NBody.days_per_year for t in time_steps],"(years)")
        pplot.plot(time_steps, KE)
        pplot.xlabel("Time " + time_unit); pplot.ylabel("Kinetic energy (J)")
        pplot.title(str(self._N) + "-body system over " + time_str)
        pplot.savefig(graphfile + "_ke.png", dpi=dpi_ke_g)


    def _write_video(self, videofile, ffmpeg_exec, time_steps, r):
        """
        Produces animation of data and saves to videofile
        """

        dpi_v = 200
        fps = 20
        margin = 1.25

        x_range = [min([min(r[i][:,0]) for i in range(self._N)]),
                   max([max(r[i][:,0]) for i in range(self._N)])]
        y_range = [min([min(r[i][:,1]) for i in range(self._N)]),
                   max([max(r[i][:,1]) for i in range(self._N)])]
        fig_vid = pplot.figure()
        ax = fig_vid.add_subplot(1,1,1)
        pplot.axis("scaled")

        def update_frame(i):
            ax.clear()
            ax.axis([a*margin for a in [x_range[0], x_range[1], y_range[0], y_range[1]]])
            for j in range(self._N):
                ax.plot(r[j][:,0][i], r[j][:,1][i],
                        "o", linewidth=0.5,
                        color=self._bodies[j].colour,
                        label=self._bodies[j].name)
            ax.legend(loc=1)
            pplot.xlabel("X (km)"); pplot.ylabel("Y (km)")
            frame_time = time_steps[i]
            time_str = str(frame_time) + " days" if frame_time < 10*NBody.days_per_year \
                        else "~" + str(int(frame_time//NBody.days_per_year)) + " years"
            pplot.title(time_str)

        vid = anim.FuncAnimation(fig_vid, update_frame, len(time_steps), interval=1000/fps)
        pplot.rcParams["animation.ffmpeg_path"] = ffmpeg_exec
        vid.save(videofile + ".mp4", anim.writers["ffmpeg"](fps=fps), dpi=dpi_v)
