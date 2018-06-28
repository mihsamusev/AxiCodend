﻿using System;
using Accord.Math;
using System.IO;
using System.Collections.Generic;

namespace AxiCodend
{
    public class AxiModelT90
    {
        //
        //
        //                  
        //            \    /\    /
        //             \  /mo\mo/
        //              \/    \/
        //               |lo   |lo   ==> current
        //              /\    /\
        //             /  \mo/  \
        //            /    \/    \
        //                     
        //
        // CLASS VARIABLES
        #region
        /*mesh characteristics*/

        private double l0;           // length of un-stretched knot       
        private double m0;           // length of un-stretched twine
        private double kl;           // stiffness of a knot
        private double km;           // stiffness of a twine

        /*codend paramters*/

        private int nx;              // amount of meshes along the codend length
        private int nr;              // amount of meshes in the codend`s circumference
        private int nc;              // moundt of meshes along the codend length contatining catch 
        private double r0;           // entrance radius

        /*operation*/

        private double P; // catch pressure
        #endregion

        /* calculation variables*/
        #region
        public readonly int dof;                    // nb of dof
        public readonly int nbc_beg = 3;	        // nb of fixed dof at entry (x1 y1 z1)
        public readonly int nbc_end = 2;            // nb of fixed dof at the end (yn zn)

        const double pi = Math.PI;
        private int ncp;                            // node number indicating start of application of the catch pressure
        private double theta, co, si;               // angle in circumference covered by 1 mesh
                                                    // sin cos of single and double theta
        private double R;                           // residue
        private double l, m;                        // length of twine and knot between 2 points, for calculation

        private double[] X0,x,y,z,r,                // initial shape, coordinates and radii vectors
                         FA,FB,FF,                  // helper force vectors (force from twines going around, backward, forward)
                         FT,FC,F;                   // twine force, catch force and total force vectors
        private double[][] dFAdX,dFBdX,dFFdX,       // helper jacobian matrices (stiffness from twines going around, backward, forward) 
                           dFTdX,dFCdX,dFdX;        // jacobian due to twine force, catch force and total force




        #endregion

        // CLASS CONSTRUCTS
        #region
        //public AxiModelT90(MeshParam _meshParam, CodendParam _codendParam, double _towSpeed)
        //{
        //    /*main inputs*/
        //    #region
        //    l0 = _meshParam.l0;
        //    m0 = _meshParam.m0;
        //    kl = _meshParam.kl;
        //    km = _meshParam.km;

        //    nx = _codendParam.nx;
        //    nr = _codendParam.nr;
        //    nc = _codendParam.nc;
        //    r0 = _codendParam.r0;

        //    P = 0.5 * 1.4 * 1025.0 * Math.Pow(_towSpeed, 2);
        //    #endregion

        //    /* sizes */
        //    dof = nx * 6 + 3;
        //    ncp = 2 * nx + 1 - 2 * nc;

        //    /*angles and trigonometric functions*/
        //    theta = 2 * pi / nr;
        //    co = Math.Cos(theta);
        //    si = Math.Sin(theta);

        //    /*force vectors and jacobian matrices*/
        //    #region
        //    X0 = new double[dof];
        //    x = new double[dof / 3];
        //    y = new double[dof / 3];
        //    z = new double[dof / 3];
        //    r = new double[dof / 3];

        //    FA = new double[dof];
        //    FB = new double[dof];
        //    FF = new double[dof];
        //    FT = new double[dof];
        //    FC = new double[dof];

        //    dFTdX = new double[dof][];
        //    dFCdX = new double[dof][];
        //    dFAdX = new double[dof][];
        //    dFBdX = new double[dof][];
        //    dFFdX = new double[dof][];
        //    dFdX = new double[dof][];
        //    for (int i = 0; i < dof; i++)
        //    {
        //        dFTdX[i] = new double[dof];
        //        dFAdX[i] = new double[dof];
        //        dFBdX[i] = new double[dof];
        //        dFFdX[i] = new double[dof];
        //        dFCdX[i] = new double[dof];
        //        dFdX[i] = new double[dof];
        //    }
        //    #endregion
        //}
        #endregion

        // CLASS METHODS

        /* for updating*/

        public void ChangeCatch(int newnc)
        {
            nc = newnc;
            ncp = 2 * nx + 1 - 2 * nc;

            Array.Clear(FC, 0, FC.Length);

            for (int i = 0; i < dof; i++)
            {
                Array.Clear(dFCdX[i], 0, dFCdX[i].Length);
            }
        }

        public void ChangeKnotLength(double newl0)
        {
            l0 = newl0;
        }

        public void ChangeSpeed(double newu)
        {
            P = 0.5 * 1.4 * 1025.0 * Math.Pow(newu, 2);
        }

        /* initial shape*/

        public double[] GetInitialShape()
        {
            X0[1] = l0 / 2;                                             // y1
            X0[2] = Math.Sqrt(Math.Pow(r0, 2) - Math.Pow((l0 / 2), 2)); // z1

            /*x for all nodes*/
            for (int i = 0; i < dof / 3; i++)
            {
                X0[3 * i] = i * (l0 + m0) / 2;
            }

            /*even nodes y and z*/
            for (int i = 1; i <= nx; i++)
            {
                X0[6 * i - 2] = r0 * Math.Sin(2 * pi / 180);          // y
                X0[6 * i - 1] = r0 * Math.Cos(2 * pi / 180);          // z
            }

            /*odd nodes y and z*/
            for (int i = 1; i <= (nx - 1); i++)
            {
                X0[6 * i + 1] = r0 * Math.Sin(2 * pi / 360);      // y
                X0[6 * i + 2] = r0 * Math.Cos(2 * pi / 360);      // z
            }

            return X0;
        }

        public double[] GetInitialShapeSmooth()
        {
            X0[1] = l0 / 2;                                             // y1
            X0[2] = Math.Sqrt(Math.Pow(r0, 2) - Math.Pow((l0 / 2), 2)); // z1

            double fStretch = 1.01;
            double L = 0.5 * pi * r0 / (fStretch * m0);
            double omega = 0;
            /*x for all nodes*/
            for (int i = 0; i < dof / 3; i++)
            {
                X0[3 * i] = i * fStretch * m0;

                if (i > dof/3 - L)
                {
                    omega += fStretch * m0 / r0;
                    X0[3 * i] = X0[3 * i - 3] + fStretch * m0 * Math.Cos(omega);
                }
            }

            omega = 0;
            
            for (int i = 1; i <= nx; i++)
            {
            /*even nodes y and z*/
                X0[6 * i - 2] = r0 * Math.Sin(2 * pi / (3 * nr));      // y
                X0[6 * i - 1] = r0 * Math.Cos(2 * pi / (3 * nr));      // z

            /*odd nodes y and z*/
                X0[6 * i + 1] = r0 * Math.Sin(2 * pi / (6 * nr));      // y
                X0[6 * i + 2] = r0 * Math.Cos(2 * pi / (6 * nr));      // z

                if(i > nx - L / 2 + 1)
                {
                    /*even nodes y and z*/
                    omega += fStretch * m0 / r0;
                    X0[6 * i - 2] = r0 * Math.Sin(2 * pi / (3 * nr)) * Math.Cos(omega);      // y
                    X0[6 * i - 1] = r0 * Math.Cos(2 * pi / (3 * nr)) * Math.Cos(omega);      // z

                    /*odd nodes y and z*/
                    omega += fStretch * m0 / r0;
                    X0[6 * i + 1] = r0 * Math.Sin(2 * pi / (6 * nr)) * Math.Cos(omega);      // y
                    X0[6 * i + 2] = r0 * Math.Cos(2 * pi / (6 * nr)) * Math.Cos(omega);      // z
                }
            }
            X0[dof - 1] = 0;
            X0[dof - 2] = 0;

            //for (int i = 0; i < dof / 3; i++)
            //{
            //    Console.WriteLine(String.Format("{0,15:N4} {1,15:N4} {2,15:N4} {3,15:N4}",
            //                                    i + 1, X0[3 * i], X0[3 * i + 1], X0[3 * i + 2]));
            //}
            return X0;
        }

        public double[] GetInitialShape(string path)
        {
            string[] lines = File.ReadAllLines(path);

            for (int i = 0; i < lines.Length; i++)
            {
                X0[i] = Convert.ToDouble(lines[i]);
            }
            return X0;
        }

        /* forces*/

        public double[] GetTwineForces(double[] X)
        {
            for (int i = 0; i < dof / 3; i++)
            {
                x[i] = X[3 * i];
                y[i] = X[3 * i + 1];
                z[i] = X[3 * i + 2];
            }

            //TENSION IN AROUND KNOTS
            #region
            /*from i-k knots - odd nodes */
            for (int i = 0; i < nx; i++)
            {
                l = 2 * y[2 * i];   // since i and k are mirrored with respect to plane XZ, length i-k is double of i coordinate

                if (l < l0)
                {
                    FA[6 * i + 1] = 0;                                  // Fy
                }
                else
                {
                    FA[6 * i + 1] = -(kl / l0) * (Math.Abs(l) - l0);    // Fy
                }
            }

            /*from i-j knots - even nodes*/
            for (int i = 1; i <= nx; i++)
            {
                l = Math.Sqrt(Math.Pow((si * z[2 * i - 1] - co * y[2 * i - 1] - y[2 * i - 1]), 2) +     // to find length i-j , j is calculated as 
                              Math.Pow((co * z[2 * i - 1] + si * y[2 * i - 1] - z[2 * i - 1]), 2));     // clockwise rotation of k by angle theta

                if (l < l0)
                {
                    FA[6 * i - 2] = 0;    // Fy
                    FA[6 * i - 1] = 0;    // Fz
                }
                else
                {
                    FA[6 * i - 2] = (kl / l0) * (l - l0) * (si * z[2 * i - 1] - co * y[2 * i - 1] - y[2 * i - 1]) / l;      // Fy
                    FA[6 * i - 1] = (kl / l0) * (l - l0) * (co * z[2 * i - 1] + si * y[2 * i - 1] - z[2 * i - 1]) / l;      // Fz
                }

            }
            #endregion

            //TENSION IN BACKWARD TWINES
            #region
            /*from twines behind even i nodes*/
            for (int i = 1; i <= nx; i++)
            {
                m = Math.Sqrt(Math.Pow((x[2 * i - 2] - x[2 * i - 1]), 2) +  // twine length behind even i node
                              Math.Pow((y[2 * i - 2] - y[2 * i - 1]), 2) +
                              Math.Pow((z[2 * i - 2] - z[2 * i - 1]), 2));
                if (m < m0)
                {
                    FB[6 * i - 3] = 0;              // Fx
                    FB[6 * i - 2] = 0;              // Fy
                    FB[6 * i - 1] = 0;              // Fz
                }
                else
                {
                    FB[6 * i - 3] = (km / m0) * (m - m0) * (x[2 * i - 2] - x[2 * i - 1]) / m;   // Fx
                    FB[6 * i - 2] = (km / m0) * (m - m0) * (y[2 * i - 2] - y[2 * i - 1]) / m;   // Fy
                    FB[6 * i - 1] = (km / m0) * (m - m0) * (z[2 * i - 2] - z[2 * i - 1]) / m;   // Fz
                }
            }

            /*from twines behind odd i nodes */
            for (int i = 1; i <= nx; i++)
            {
                m = Math.Sqrt(Math.Pow((x[2 * i - 1] - x[2 * i]), 2) +  // twine length behind odd i node
                              Math.Pow((y[2 * i - 1] - y[2 * i]), 2) +
                              Math.Pow((z[2 * i - 1] - z[2 * i]), 2));
                if (m < m0)
                {
                    FB[6 * i + 0] = 0;  // Fx
                    FB[6 * i + 1] = 0;  // Fy
                    FB[6 * i + 2] = 0;  // Fz
                }
                else
                {
                    FB[6 * i + 0] = (km / m0) * (m - m0) * (x[2 * i - 1] - x[2 * i]) / m;   // Fx
                    FB[6 * i + 1] = (km / m0) * (m - m0) * (y[2 * i - 1] - y[2 * i]) / m;   // Fy
                    FB[6 * i + 2] = (km / m0) * (m - m0) * (z[2 * i - 1] - z[2 * i]) / m;   // Fz
                }
            }
            #endregion

            //TENSION IN FORWARD TWINES
            #region
            /*from twines after even i nodes*/
            for (int i = 1; i <= nx; i++)
            {
                m = Math.Sqrt(Math.Pow((x[2 * i] - x[2 * i - 1]), 2) +  // twine length after even i node
                              Math.Pow((y[2 * i] - y[2 * i - 1]), 2) +
                              Math.Pow((z[2 * i] - z[2 * i - 1]), 2));
                if (m < m0)
                {
                    FF[6 * i - 3] = 0;  // Fx
                    FF[6 * i - 2] = 0;  // Fy
                    FF[6 * i - 1] = 0;  // Fz
                }
                else
                {
                    FF[6 * i - 3] = (km / m0) * (m - m0) * (x[2 * i] - x[2 * i - 1]) / m;  // Fx
                    FF[6 * i - 2] = (km / m0) * (m - m0) * (y[2 * i] - y[2 * i - 1]) / m;  // Fy
                    FF[6 * i - 1] = (km / m0) * (m - m0) * (z[2 * i] - z[2 * i - 1]) / m;  // Fz
                }
            }

            /*from twines after odd i nodes */
            for (int i = 0; i <= (nx - 1); i++)
            {
                m = Math.Sqrt(Math.Pow(x[2 * i + 1] - x[2 * i], 2) +  // twine length after odd i node
                              Math.Pow(y[2 * i + 1] - y[2 * i], 2) +
                              Math.Pow(z[2 * i + 1] - z[2 * i], 2));
                if (m < m0)
                {
                    FF[6 * i + 0] = 0;  // Fx
                    FF[6 * i + 1] = 0;  // Fy
                    FF[6 * i + 2] = 0;  // Fz
                }
                else
                {
                    FF[6 * i + 0] = (km / m0) * (m - m0) * (x[2 * i + 1] - x[2 * i]) / m;  // Fx
                    FF[6 * i + 1] = (km / m0) * (m - m0) * (y[2 * i + 1] - y[2 * i]) / m;  // Fy
                    FF[6 * i + 2] = (km / m0) * (m - m0) * (z[2 * i + 1] - z[2 * i]) / m;  // Fz
                }
            }
            #endregion

            FT = FA.Add(FF.Add(FB));    // FT = FA + FF + FB collect all contributions at each node

            return FT;
        }

        public double[] GetCatchForces(double[] X)
        {
            Array.Clear(FF, 0, FF.Length);
            Array.Clear(FB, 0, FB.Length);

            for (int i = 0; i < dof / 3; i++)
            {
                x[i] = X[3 * i];
                y[i] = X[3 * i + 1];
                z[i] = X[3 * i + 2];
                r[i] = Math.Sqrt(Math.Pow(y[i], 2) + Math.Pow(z[i], 2));
            }

            //TENSION DUE TO CATCH PRESSURE BACKWARD
            #region
            for (int i = ncp + 1; i <= 2 * nx + 1; i++)
            {                            
                if (i == 2 * nx + 1)   // very last point
                {
                    FB[3 * i - 3] = pi * 0.25 * P / nr * (Math.Pow((r[i - 2]), 2) - Math.Pow((r[i - 1]), 2));
                    FB[3 * i - 2] = pi * 0.25 * P / nr * (r[i - 2] + r[i - 1]) * (x[i - 1] - x[i - 2]) / Math.Pow(2, 0.5);
                    FB[3 * i - 1] = pi * 0.25 * P / nr * (r[i - 2] + r[i - 1]) * (x[i - 1] - x[i - 2]) / Math.Pow(2, 0.5);
                }
                else                     // pressure backward starting from npp+1
                {
                FB[3 * i - 3] = pi * 0.25 * P / nr * (Math.Pow((r[i - 2]), 2) - Math.Pow((r[i - 1]), 2));
                FB[3 * i - 2] = pi * 0.25 * P / nr * (r[i - 2] + r[i - 1]) * y[i - 1] / r[i - 1] * (x[i - 1] - x[i - 2]);
                FB[3 * i - 1] = pi * 0.25 * P / nr * (r[i - 2] + r[i - 1]) * z[i - 1] / r[i - 1] * (x[i - 1] - x[i - 2]);
                }
            }
            #endregion

            //TENSION DUE TO CATCH PRESSURE FORWARD
            #region
            for (int i = ncp; i <= 2 * nx; i++)
            {
                FF[3 * i - 3] = pi * 0.25 * P / nr * (Math.Pow((r[i - 1]), 2) - Math.Pow((r[i]), 2));
                FF[3 * i - 2] = pi * 0.25 * P / nr * (r[i - 1] + r[i]) * y[i - 1] / r[i - 1] * (x[i] - x[i - 1]);
                FF[3 * i - 1] = pi * 0.25 * P / nr * (r[i - 1] + r[i]) * z[i - 1] / r[i - 1] * (x[i] - x[i - 1]);
            }
            #endregion

            FC = FF.Add(FB); // FC = FF + FB collect all contributions at each node

            return FC;
        }

        public double[] GetTotalForces(double[] X)
        {
            FT = GetTwineForces(X);
            FC = GetCatchForces(X);
         
            F = FT.Add(FC);             /*Add elastic twine forces and catch presure forces*/

            return F;
        }

        public double ForceResidue(double[] F)
        {
            R = 0;

            for (int i = nbc_beg; i < dof - nbc_end; i++)
            {
                R += Math.Pow(F[i], 2);                      // residue as sum of squares
            }
            R = Math.Sqrt(R);
            return R;
        }

        /* jacobian and stiffness*/

        private void BandAdd2(ref double[][] AB,double[][] A, double[][] B)
        {
            int shift = 0;                    // Faster way is to add only the band terms ignoring zeros
            int width = 0;

            for (int i = 0; i < dof; i++)
            {
                if (i < 3)
                {
                    width = 6;
                    shift = 0;
                }
                else if (i > dof - 4)
                {
                    width = 6;
                    shift = dof - width;
                }
                else
                {
                    shift = shift = 3 * (i / 3) - 3; // integer division (no remainder)
                    width = 9;
                }

                for (int j = 0; j < width; j++)
                {
                    AB[i][j + shift] =   A[i][j + shift]
                                       + B[i][j + shift];
                }
            }
        }

        private void BandAdd3(ref double[][] ABC, double[][] A, double[][] B, double[][] C)
        {
            int shift = 0;                    // Faster way is to add only the band terms ignoring zeros
            int width = 0;

            for (int i = 0; i < dof; i++)
            {
                if (i < 3)
                {
                    width = 6;
                    shift = 0;
                }
                else if (i > dof - 4)
                {
                    width = 6;
                    shift = dof - width;
                }
                else
                {
                    shift = shift = 3 * (i / 3) - 3; // integer division (no remainder)
                    width = 9;
                }

                for (int j = 0; j < width; j++)
                {
                    ABC[i][j + shift] =  A[i][j + shift]
                                       + B[i][j + shift]
                                       + C[i][j + shift];
                }
            }
        }

        public double[][] GetTwineJacobian(double[] X)
        {
            for (int i = 0; i < dof / 3; i++)
            {
                x[i] = X[3 * i];
                y[i] = X[3 * i + 1];
                z[i] = X[3 * i + 2];
            }

            //TENSION IN AROUND KNOTS
            #region
            /*from i-k knots - odd nodes */
            for (int i = 0; i < nx; i++)
            {
                l = 2 * y[2 * i];   // since i and k are mirrored with respect to plane XZ, length i-k is double of i coordinate

                if (l < l0)
                {
                    dFAdX[6 * i + 1][6 * i + 1] = 0;                // Fy with respect to y
                }
                else
                {
                    dFAdX[6 * i + 1][6 * i + 1] = -2 * kl / l0;     // Fy with respect to y
                }
            }

            /*from i-j knots - even nodes*/
            for (int i = 1; i <= nx; i++)
            {
                l = Math.Sqrt(Math.Pow((si * z[2 * i - 1] - co * y[2 * i - 1] - y[2 * i - 1]), 2) +     // to find length i-j , j is calculated as 
                              Math.Pow((co * z[2 * i - 1] + si * y[2 * i - 1] - z[2 * i - 1]), 2));     // clockwise rotation of k by angle theta

                if (l < l0)
                {
                    // Fy with respect to y and z
                    dFAdX[6 * i - 2][6 * i - 2] = 0;    
                    dFAdX[6 * i - 2][6 * i - 1] = 0;
                    // Fz with respect to y and z
                    dFAdX[6 * i - 1][6 * i - 2] = 0;
                    dFAdX[6 * i - 1][6 * i - 1] = 0;
                }
                else
                {
                    // Fy with respect to y and z
                    dFAdX[6 * i - 2][6 * i - 2] = ((-co - 1) * kl * (l - l0)) / (l0 * l) - (0.5 * kl * (-co * y[2 * i - 1] - y[2 * i - 1] + si * z[2 * i - 1]) * (2 * si * (si * y[2 * i - 1] + co * z[2 * i - 1] - z[2 * i - 1]) + 2 * (-co - 1) * (-co * y[2 * i - 1] - y[2 * i - 1] + si * z[2 * i - 1])) * (l - l0)) / (l0 * Math.Pow(l, 3)) + (0.5 * kl * (-co * y[2 * i - 1] - y[2 * i - 1] + si * z[2 * i - 1]) * (2 * si * (si * y[2 * i - 1] + co * z[2 * i - 1] - z[2 * i - 1]) + 2 * (-co - 1) * (-co * y[2 * i - 1] - y[2 * i - 1] + si * z[2 * i - 1]))) / (l0 * Math.Pow(l, 2));
                    dFAdX[6 * i - 2][6 * i - 1] = (si * kl * (l - l0)) / (l0 * l) - (0.5 * kl * (-co * y[2 * i - 1] - y[2 * i - 1] + si * z[2 * i - 1]) * (2 * (co - 1) * (si * y[2 * i - 1] + co * z[2 * i - 1] - z[2 * i - 1]) + 2 * si * (-co * y[2 * i - 1] - y[2 * i - 1] + si * z[2 * i - 1])) * (l - l0)) / (l0 * Math.Pow(l, 3)) + (0.5 * kl * (-co * y[2 * i - 1] - y[2 * i - 1] + si * z[2 * i - 1]) * (2 * (co - 1) * (si * y[2 * i - 1] + co * z[2 * i - 1] - z[2 * i - 1]) + 2 * si * (-co * y[2 * i - 1] - y[2 * i - 1] + si * z[2 * i - 1]))) / (l0 * Math.Pow(l, 2));
                    // Fz with respect to y and z
                    dFAdX[6 * i - 1][6 * i - 2] = (si * kl * (l - l0)) / (l0 * l) - (0.5 * kl * (si * y[2 * i - 1] + co * z[2 * i - 1] - z[2 * i - 1]) * (2 * si * (si * y[2 * i - 1] + co * z[2 * i - 1] - z[2 * i - 1]) + 2 * (-co - 1) * (-co * y[2 * i - 1] - y[2 * i - 1] + si * z[2 * i - 1])) * (l - l0)) / (l0 * Math.Pow(l, 3)) + (0.5 * kl * (si * y[2 * i - 1] + co * z[2 * i - 1] - z[2 * i - 1]) * (2 * si * (si * y[2 * i - 1] + co * z[2 * i - 1] - z[2 * i - 1]) + 2 * (-co - 1) * (-co * y[2 * i - 1] - y[2 * i - 1] + si * z[2 * i - 1]))) / (l0 * Math.Pow(l, 2));
                    dFAdX[6 * i - 1][6 * i - 1] = ((co - 1) * kl * (l - l0)) / (l0 * l) - (0.5 * kl * (si * y[2 * i - 1] + co * z[2 * i - 1] - z[2 * i - 1]) * (2 * (co - 1) * (si * y[2 * i - 1] + co * z[2 * i - 1] - z[2 * i - 1]) + 2 * si * (-co * y[2 * i - 1] - y[2 * i - 1] + si * z[2 * i - 1])) * (l - l0)) / (l0 * Math.Pow(l, 3)) + (0.5 * kl * (si * y[2 * i - 1] + co * z[2 * i - 1] - z[2 * i - 1]) * (2 * (co - 1) * (si * y[2 * i - 1] + co * z[2 * i - 1] - z[2 * i - 1]) + 2 * si * (-co * y[2 * i - 1] - y[2 * i - 1] + si * z[2 * i - 1]))) / (l0 * Math.Pow(l, 2));
                }
            }
            #endregion

            //TENSION IN BACKWARD TWINES
            #region
            /*from twines behind even i nodes*/
            for (int i = 1; i <= nx; i++)
            {
                m = Math.Sqrt(Math.Pow((x[2 * i - 2] - x[2 * i - 1]), 2) +  // twine length behind even i node
                              Math.Pow((y[2 * i - 2] - y[2 * i - 1]), 2) +
                              Math.Pow((z[2 * i - 2] - z[2 * i - 1]), 2));
                if (m < m0)
                {
                    for (int j = 6; j > 0; j--)
                    {
                        // Fx with respect to x y z (this and previous node)
                        dFBdX[6 * i - 3][6 * i - j] = 0;
                        // Fy with respect to x y z (this and previous node)
                        dFBdX[6 * i - 2][6 * i - j] = 0;
                        // Fz with respect to x y z (this and previous node)
                        dFBdX[6 * i - 1][6 * i - j] = 0;
                    }
                }
                else
                {
                    // Fx with respect to x y z (this and previous node)
                    dFBdX[6 * i - 3][6 * i - 6] = -(-km * (m - m0) / m0 / m + km * Math.Pow((x[2 * i - 1] - x[2 * i - 2]), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(((x[2 * i - 1] - x[2 * i - 2])), 2) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i - 3][6 * i - 5] = -(km * (x[2 * i - 1] - x[2 * i - 2]) * (y[2 * i - 1] - y[2 * i - 2]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i - 1] - x[2 * i - 2]) * (y[2 * i - 1] - y[2 * i - 2]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i - 3][6 * i - 4] = -(km * (x[2 * i - 1] - x[2 * i - 2]) * (z[2 * i - 1] - z[2 * i - 2]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i - 1] - x[2 * i - 2]) * (z[2 * i - 1] - z[2 * i - 2]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i - 3][6 * i - 3] = (-km * (m - m0) / m0 / m + km * Math.Pow(((x[2 * i - 1] - x[2 * i - 2])), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(((x[2 * i - 1] - x[2 * i - 2])), 2) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i - 3][6 * i - 2] = (km * (x[2 * i - 1] - x[2 * i - 2]) * (y[2 * i - 1] - y[2 * i - 2]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i - 1] - x[2 * i - 2]) * (y[2 * i - 1] - y[2 * i - 2]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i - 3][6 * i - 1] = -(-(km * (x[2 * i - 1] - x[2 * i - 2]) * (z[2 * i - 1] - z[2 * i - 2]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i - 1] - x[2 * i - 2]) * (z[2 * i - 1] - z[2 * i - 2]) / m0 / Math.Pow(m, 2)));
                    // Fy with respect to x y z (this and previous node)
                    dFBdX[6 * i - 2][6 * i - 6] = -(-(km * (x[2 * i - 2] - x[2 * i - 1]) * (y[2 * i - 1] - y[2 * i - 2]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i - 2] - x[2 * i - 1]) * (y[2 * i - 1] - y[2 * i - 2]) / m0 / Math.Pow(m, 2)));
                    dFBdX[6 * i - 2][6 * i - 5] = -(-km * (m - m0) / m0 / m + km * Math.Pow(((y[2 * i - 1] - y[2 * i - 2])), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(((y[2 * i - 1] - y[2 * i - 2])), 2) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i - 2][6 * i - 4] = -(km * ((y[2 * i - 1] - y[2 * i - 2])) * (z[2 * i - 1] - z[2 * i - 2]) * (m - m0) / m0 / Math.Pow(m, 3) - km * ((y[2 * i - 1] - y[2 * i - 2])) * (z[2 * i - 1] - z[2 * i - 2]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i - 2][6 * i - 3] = -(km * (x[2 * i - 2] - x[2 * i - 1]) * (y[2 * i - 1] - y[2 * i - 2]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i - 2] - x[2 * i - 1]) * (y[2 * i - 1] - y[2 * i - 2]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i - 2][6 * i - 2] = (-km * (m - m0) / m0 / m + km * Math.Pow(((y[2 * i - 1] - y[2 * i - 2])), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(((y[2 * i - 1] - y[2 * i - 2])), 2) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i - 2][6 * i - 1] = (km * ((y[2 * i - 1] - y[2 * i - 2])) * (z[2 * i - 1] - z[2 * i - 2]) * (m - m0) / m0 / Math.Pow(m, 3) - km * ((y[2 * i - 1] - y[2 * i - 2])) * (z[2 * i - 1] - z[2 * i - 2]) / m0 / Math.Pow(m, 2));
                    // Fz with respect to x y z (this and previous node)
                    dFBdX[6 * i - 1][6 * i - 6] = (km * (x[2 * i - 2] - x[2 * i - 1]) * (z[2 * i - 1] - z[2 * i - 2]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i - 2] - x[2 * i - 1]) * (z[2 * i - 1] - z[2 * i - 2]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i - 1][6 * i - 5] = -(km * ((y[2 * i - 1] - y[2 * i - 2])) * (z[2 * i - 1] - z[2 * i - 2]) * (m - m0) / m0 / Math.Pow(m, 3) - km * ((y[2 * i - 1] - y[2 * i - 2])) * (z[2 * i - 1] - z[2 * i - 2]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i - 1][6 * i - 4] = -(-km * (m - m0) / m0 / m + km * Math.Pow(((z[2 * i - 1] - z[2 * i - 2])), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(((z[2 * i - 1] - z[2 * i - 2])), 2) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i - 1][6 * i - 3] = -(km * (x[2 * i - 2] - x[2 * i - 1]) * (z[2 * i - 1] - z[2 * i - 2]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i - 2] - x[2 * i - 1]) * (z[2 * i - 1] - z[2 * i - 2]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i - 1][6 * i - 2] = (km * ((y[2 * i - 1] - y[2 * i - 2])) * (z[2 * i - 1] - z[2 * i - 2]) * (m - m0) / m0 / Math.Pow(m, 3) - km * ((y[2 * i - 1] - y[2 * i - 2])) * (z[2 * i - 1] - z[2 * i - 2]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i - 1][6 * i - 1] = (-km * (m - m0) / m0 / m + km * Math.Pow(((z[2 * i - 1] - z[2 * i - 2])), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(((z[2 * i - 1] - z[2 * i - 2])), 2) / m0 / Math.Pow(m, 2));
                }
            }

            /*from twines behind odd i nodes */
            for (int i = 1; i <= nx; i++)
            {
                m = Math.Sqrt(Math.Pow((x[2 * i - 1] - x[2 * i]), 2) +  // twine length behind odd i node
                              Math.Pow((y[2 * i - 1] - y[2 * i]), 2) +
                              Math.Pow((z[2 * i - 1] - z[2 * i]), 2));
                if (m < m0)
                {
                    for (int j = 3; j > -3; j--)
                    {
                        // Fx with respect to x y z (this and previous node)
                        dFBdX[6 * i + 0][6 * i - j] = 0;
                        // Fy with respect to x y z (this and previous node)
                        dFBdX[6 * i + 1][6 * i - j] = 0;
                        // Fz with respect to x y z (this and previous node)
                        dFBdX[6 * i + 2][6 * i - j] = 0;
                    }
                }
                else
                {
                    // Fx with respect to x y z (this and previous node)
                    dFBdX[6 * i + 0][6 * i - 3] = -(-km * (m - m0) / m0 / m + km * Math.Pow((x[2 * i] - x[2 * i - 1]), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow((x[2 * i] - x[2 * i - 1]), 2) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i + 0][6 * i - 2] = -(km * (x[2 * i] - x[2 * i - 1]) * (y[2 * i] - y[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i] - x[2 * i - 1]) * (y[2 * i] - y[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i + 0][6 * i - 1] = -(km * (x[2 * i] - x[2 * i - 1]) * (z[2 * i] - z[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i] - x[2 * i - 1]) * (z[2 * i] - z[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i + 0][6 * i - 0] = (-km * (m - m0) / m0 / m + km * Math.Pow((x[2 * i] - x[2 * i - 1]), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow((x[2 * i] - x[2 * i - 1]), 2) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i + 0][6 * i + 1] = (km * (x[2 * i] - x[2 * i - 1]) * (y[2 * i] - y[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i] - x[2 * i - 1]) * (y[2 * i] - y[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i + 0][6 * i + 2] = (km * (x[2 * i] - x[2 * i - 1]) * (z[2 * i] - z[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i] - x[2 * i - 1]) * (z[2 * i] - z[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    // Fy with respect to x y z (this and previous node)
                    dFBdX[6 * i + 1][6 * i - 3] = -(km * (x[2 * i] - x[2 * i - 1]) * (y[2 * i] - y[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i] - x[2 * i - 1]) * (y[2 * i] - y[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i + 1][6 * i - 2] = -(-km * (m - m0) / m0 / m + km * Math.Pow(((y[2 * i] - y[2 * i - 1])), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(((y[2 * i] - y[2 * i - 1])), 2) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i + 1][6 * i - 1] = -(km * ((y[2 * i] - y[2 * i - 1])) * (z[2 * i] - z[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * ((y[2 * i] - y[2 * i - 1])) * (z[2 * i] - z[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i + 1][6 * i - 0] = (km * (x[2 * i] - x[2 * i - 1]) * (y[2 * i] - y[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i] - x[2 * i - 1]) * (y[2 * i] - y[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i + 1][6 * i + 1] = (-km * (m - m0) / m0 / m + km * Math.Pow(((y[2 * i] - y[2 * i - 1])), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(((y[2 * i] - y[2 * i - 1])), 2) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i + 1][6 * i + 2] = (km * ((y[2 * i] - y[2 * i - 1])) * (z[2 * i] - z[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * ((y[2 * i] - y[2 * i - 1])) * (z[2 * i] - z[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    // Fz with respect to x y z (this and previous node)
                    dFBdX[6 * i + 2][6 * i - 3] = -(km * (x[2 * i] - x[2 * i - 1]) * (z[2 * i] - z[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i] - x[2 * i - 1]) * (z[2 * i] - z[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i + 2][6 * i - 2] = -(km * ((y[2 * i] - y[2 * i - 1])) * (z[2 * i] - z[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * ((y[2 * i] - y[2 * i - 1])) * (z[2 * i] - z[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i + 2][6 * i - 1] = -(-km * (m - m0) / m0 / m + km * Math.Pow(((z[2 * i] - z[2 * i - 1])), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(((z[2 * i] - z[2 * i - 1])), 2) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i + 2][6 * i - 0] = (km * (x[2 * i] - x[2 * i - 1]) * (z[2 * i] - z[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i] - x[2 * i - 1]) * (z[2 * i] - z[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i + 2][6 * i + 1] = (km * ((y[2 * i] - y[2 * i - 1])) * (z[2 * i] - z[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * ((y[2 * i] - y[2 * i - 1])) * (z[2 * i] - z[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFBdX[6 * i + 2][6 * i + 2] = (-km * (m - m0) / m0 / m + km * Math.Pow(((z[2 * i] - z[2 * i - 1])), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(((z[2 * i] - z[2 * i - 1])), 2) / m0 / Math.Pow(m, 2));
                }
            }
            #endregion

            //TENSION IN FORWARD TWINES
            #region
            /*from twines after even i nodes*/
            for (int i = 1; i <= nx; i++)
            {
                m = Math.Sqrt(Math.Pow((x[2 * i] - x[2 * i - 1]), 2) +  // twine length after even i node
                              Math.Pow((y[2 * i] - y[2 * i - 1]), 2) +
                              Math.Pow((z[2 * i] - z[2 * i - 1]), 2));
                if (m < m0)
                {
                    for (int j = 3; j > -3; j--)
                    {
                        // Fx with respect to x y z (this and next node)
                        dFFdX[6 * i - 3][6 * i - j] = 0;
                        // Fy with respect to x y z (this and next node)
                        dFFdX[6 * i - 2][6 * i - j] = 0;
                        // Fz with respect to x y z (this and next node)
                        dFFdX[6 * i - 1][6 * i - j] = 0;
                    }
                }
                else
                {
                    // Fx with respect to x y z (this and next node)
                    dFFdX[6 * i - 3][6 * i - 3] = (-km * (m - m0) / m0 / m + km * Math.Pow((x[2 * i] - x[2 * i - 1]), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow((x[2 * i] - x[2 * i - 1]), 2) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i - 3][6 * i - 2] = (km * (x[2 * i] - x[2 * i - 1]) * (y[2 * i] - y[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i] - x[2 * i - 1]) * (y[2 * i] - y[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i - 3][6 * i - 1] = (km * (x[2 * i] - x[2 * i - 1]) * (z[2 * i] - z[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i] - x[2 * i - 1]) * (z[2 * i] - z[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i - 3][6 * i - 0] = -(-km * (m - m0) / m0 / m + km * Math.Pow((x[2 * i] - x[2 * i - 1]), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow((x[2 * i] - x[2 * i - 1]), 2) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i - 3][6 * i + 1] = -(km * (x[2 * i] - x[2 * i - 1]) * (y[2 * i] - y[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i] - x[2 * i - 1]) * (y[2 * i] - y[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i - 3][6 * i + 2] = -(km * (x[2 * i] - x[2 * i - 1]) * (z[2 * i] - z[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i] - x[2 * i - 1]) * (z[2 * i] - z[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    // Fy with respect to x y z (this and next node)
                    dFFdX[6 * i - 2][6 * i - 3] = -(km * (x[2 * i - 1] - x[2 * i]) * (y[2 * i] - y[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i - 1] - x[2 * i]) * (y[2 * i] - y[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i - 2][6 * i - 2] = (-km * (m - m0) / m0 / m + km * Math.Pow(((y[2 * i] - y[2 * i - 1])), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(((y[2 * i] - y[2 * i - 1])), 2) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i - 2][6 * i - 1] = (km * ((y[2 * i] - y[2 * i - 1])) * (z[2 * i] - z[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * ((y[2 * i] - y[2 * i - 1])) * (z[2 * i] - z[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i - 2][6 * i - 0] = (km * (x[2 * i - 1] - x[2 * i]) * (y[2 * i] - y[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i - 1] - x[2 * i]) * (y[2 * i] - y[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i - 2][6 * i + 1] = -(-km * (m - m0) / m0 / m + km * Math.Pow(((y[2 * i] - y[2 * i - 1])), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(((y[2 * i] - y[2 * i - 1])), 2) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i - 2][6 * i + 2] = -(km * ((y[2 * i] - y[2 * i - 1])) * (z[2 * i] - z[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * ((y[2 * i] - y[2 * i - 1])) * (z[2 * i] - z[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    // Fz with respect to x y z (this and next node)
                    dFFdX[6 * i - 1][6 * i - 3] = -(km * (x[2 * i - 1] - x[2 * i]) * (z[2 * i] - z[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i - 1] - x[2 * i]) * (z[2 * i] - z[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i - 1][6 * i - 2] = (km * ((y[2 * i] - y[2 * i - 1])) * (z[2 * i] - z[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * ((y[2 * i] - y[2 * i - 1])) * (z[2 * i] - z[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i - 1][6 * i - 1] = (-km * (m - m0) / m0 / m + km * Math.Pow(((z[2 * i] - z[2 * i - 1])), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(((z[2 * i] - z[2 * i - 1])), 2) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i - 1][6 * i - 0] = (km * (x[2 * i - 1] - x[2 * i]) * (z[2 * i] - z[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i - 1] - x[2 * i]) * (z[2 * i] - z[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i - 1][6 * i + 1] = -(km * ((y[2 * i] - y[2 * i - 1])) * (z[2 * i] - z[2 * i - 1]) * (m - m0) / m0 / Math.Pow(m, 3) - km * ((y[2 * i] - y[2 * i - 1])) * (z[2 * i] - z[2 * i - 1]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i - 1][6 * i + 2] = -(-km * (m - m0) / m0 / m + km * Math.Pow(((z[2 * i] - z[2 * i - 1])), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(((z[2 * i] - z[2 * i - 1])), 2) / m0 / Math.Pow(m, 2));
                }
            }

            /*from twines after odd i nodes */
            for (int i = 0; i <= (nx - 1); i++)
            {
                m = Math.Sqrt(Math.Pow(x[2 * i + 1] - x[2 * i], 2) +  // twine length after odd i node
                              Math.Pow(y[2 * i + 1] - y[2 * i], 2) +
                              Math.Pow(z[2 * i + 1] - z[2 * i], 2));
                if (m < m0)
                {
                    for (int j = 0; j < 6; j++)
                    {
                        // Fx with respect to x y z (this and next node)
                        dFFdX[6 * i + 0][6 * i + j] = 0;
                        // Fy with respect to x y z (this and next node)
                        dFFdX[6 * i + 1][6 * i + j] = 0;
                        // Fz with respect to x y z (this and next node)
                        dFFdX[6 * i + 2][6 * i + j] = 0;
                    }
                }
                else
                {
                    // Fx with respect to x y z (this and next node)
                    dFFdX[6 * i + 0][6 * i + 0] = (-km * (m - m0) / m0 / m + km * Math.Pow(x[2 * i + 1] - x[2 * i], 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * (Math.Pow(x[2 * i + 1] - x[2 * i], 2) / m0 / Math.Pow(m, 2)));
                    dFFdX[6 * i + 0][6 * i + 1] = (km * (x[2 * i + 1] - x[2 * i]) * (y[2 * i + 1] - y[2 * i]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i + 1] - x[2 * i]) * (y[2 * i + 1] - y[2 * i]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i + 0][6 * i + 2] = (km * (x[2 * i + 1] - x[2 * i]) * (z[2 * i + 1] - z[2 * i]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i + 1] - x[2 * i]) * (z[2 * i + 1] - z[2 * i]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i + 0][6 * i + 3] = -(-km * (m - m0) / m0 / m + km * Math.Pow(x[2 * i + 1] - x[2 * i], 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(x[2 * i + 1] - x[2 * i], 2) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i + 0][6 * i + 4] = -(km * (x[2 * i + 1] - x[2 * i]) * (y[2 * i + 1] - y[2 * i]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i + 1] - x[2 * i]) * (y[2 * i + 1] - y[2 * i]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i + 0][6 * i + 5] = -(km * (x[2 * i + 1] - x[2 * i]) * (z[2 * i + 1] - z[2 * i]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i + 1] - x[2 * i]) * (z[2 * i + 1] - z[2 * i]) / m0 / Math.Pow(m, 2));
                    // Fy with respect to x y z (this and next node)
                    dFFdX[6 * i + 1][6 * i + 0] = (km * (x[2 * i + 1] - x[2 * i]) * (y[2 * i + 1] - y[2 * i]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i + 1] - x[2 * i]) * (y[2 * i + 1] - y[2 * i]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i + 1][6 * i + 1] = (-km * (m - m0) / m0 / m + km * Math.Pow(((y[2 * i + 1] - y[2 * i])), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(((y[2 * i + 1] - y[2 * i])), 2) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i + 1][6 * i + 2] = (km * ((y[2 * i + 1] - y[2 * i])) * (z[2 * i + 1] - z[2 * i]) * (m - m0) / m0 / Math.Pow(m, 3) - km * ((y[2 * i + 1] - y[2 * i])) * (z[2 * i + 1] - z[2 * i]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i + 1][6 * i + 3] = -(km * (x[2 * i + 1] - x[2 * i]) * (y[2 * i + 1] - y[2 * i]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i + 1] - x[2 * i]) * (y[2 * i + 1] - y[2 * i]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i + 1][6 * i + 4] = -(-km * (m - m0) / m0 / m + km * Math.Pow(((y[2 * i + 1] - y[2 * i])), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(((y[2 * i + 1] - y[2 * i])), 2) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i + 1][6 * i + 5] = -(km * ((y[2 * i + 1] - y[2 * i])) * (z[2 * i + 1] - z[2 * i]) * (m - m0) / m0 / Math.Pow(m, 3) - km * ((y[2 * i + 1] - y[2 * i])) * (z[2 * i + 1] - z[2 * i]) / m0 / Math.Pow(m, 2));
                    // Fz with respect to x y z (this and next node)
                    dFFdX[6 * i + 2][6 * i + 0] = (km * (x[2 * i + 1] - x[2 * i]) * (z[2 * i + 1] - z[2 * i]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i + 1] - x[2 * i]) * (z[2 * i + 1] - z[2 * i]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i + 2][6 * i + 1] = (km * ((y[2 * i + 1] - y[2 * i])) * (z[2 * i + 1] - z[2 * i]) * (m - m0) / m0 / Math.Pow(m, 3) - km * ((y[2 * i + 1] - y[2 * i])) * (z[2 * i + 1] - z[2 * i]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i + 2][6 * i + 2] = (-km * (m - m0) / m0 / m + km * Math.Pow(((z[2 * i + 1] - z[2 * i])), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(((z[2 * i + 1] - z[2 * i])), 2) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i + 2][6 * i + 3] = -(km * (x[2 * i + 1] - x[2 * i]) * (z[2 * i + 1] - z[2 * i]) * (m - m0) / m0 / Math.Pow(m, 3) - km * (x[2 * i + 1] - x[2 * i]) * (z[2 * i + 1] - z[2 * i]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i + 2][6 * i + 4] = -(km * ((y[2 * i + 1] - y[2 * i])) * (z[2 * i + 1] - z[2 * i]) * (m - m0) / m0 / Math.Pow(m, 3) - km * ((y[2 * i + 1] - y[2 * i])) * (z[2 * i + 1] - z[2 * i]) / m0 / Math.Pow(m, 2));
                    dFFdX[6 * i + 2][6 * i + 5] = -(-km * (m - m0) / m0 / m + km * Math.Pow(((z[2 * i + 1] - z[2 * i])), 2) * (m - m0) / m0 / Math.Pow(m, 3) - km * Math.Pow(((z[2 * i + 1] - z[2 * i])), 2) / m0 / Math.Pow(m, 2));
                }
            }
            #endregion
            
            //dFTdX = dFAdX.Add(dFFdX.Add(dFBdX));    // dFTdX = dFAdX + dFFdX + dFBdX collect all contributions at each node
            BandAdd3(ref dFTdX, dFAdX, dFFdX, dFBdX);
            return dFTdX;
        }

        public double[][] GetCatchJacobian(double[] X)
        {
            for (int i = 0; i < dof; i++)
            {
                Array.Clear(dFFdX[i], 0, dFFdX[i].Length);
                Array.Clear(dFBdX[i], 0, dFBdX[i].Length);
            }

            for (int i = 0; i < dof / 3; i++)
            {
                x[i] = X[3 * i];
                y[i] = X[3 * i + 1];
                z[i] = X[3 * i + 2];
                r[i] = Math.Sqrt(Math.Pow(y[i], 2) + Math.Pow(z[i], 2));
            }

            //TENSION DUE TO CATCH PRESSURE BACKWARD
            #region
            for (int i = ncp + 1; i <= 2 * nx + 1; i++)
            {
                if (i == 2 * nx + 1)   // very last point
                {
                    // Fx with respect to x y z (this and previous node)
                    dFBdX[3 * i - 3][3 * i - 6] = 0;
                    dFBdX[3 * i - 3][3 * i - 5] = pi * 0.25 * P / nr * 2 * y[i - 2];
                    dFBdX[3 * i - 3][3 * i - 4] = pi * 0.25 * P / nr * 2 * z[i - 2];
                    dFBdX[3 * i - 3][3 * i - 3] = 0;
                    dFBdX[3 * i - 3][3 * i - 2] = pi * 0.25 * P / nr * (-2 * y[i - 1]);
                    dFBdX[3 * i - 3][3 * i - 1] = pi * 0.25 * P / nr * (-2 * z[i - 1]);
                    // Fy with respect to x y z (this and previous node)
                    dFBdX[3 * i - 2][3 * i - 6] = 1 / Math.Pow(2, 0.5) * pi * 0.25 * P / nr * (-(r[i - 2] + r[i - 1]));
                    dFBdX[3 * i - 2][3 * i - 5] = 1 / Math.Pow(2, 0.5) * pi * 0.25 * P / nr * ((x[i - 1] - x[i - 2]) * y[i - 2] / r[i - 2]);
                    dFBdX[3 * i - 2][3 * i - 4] = 1 / Math.Pow(2, 0.5) * pi * 0.25 * P / nr * ((x[i - 1] - x[i - 2]) * z[i - 2] / r[i - 2]);
                    dFBdX[3 * i - 2][3 * i - 3] = 1 / Math.Pow(2, 0.5) * pi * 0.25 * P / nr * (r[i - 1] + r[i - 2]);
                    dFBdX[3 * i - 2][3 * i - 2] = 0;
                    dFBdX[3 * i - 2][3 * i - 1] = 0;
                    // Fz with respect to x y z (this and previous node)
                    dFBdX[3 * i - 1][3 * i - 6] = 1 / Math.Pow(2, 0.5) * pi * 0.25 * P / nr * (-(r[i - 2] + r[i - 1]));
                    dFBdX[3 * i - 1][3 * i - 5] = 1 / Math.Pow(2, 0.5) * pi * 0.25 * P / nr * ((x[i - 1] - x[i - 2]) * y[i - 2] / r[i - 2]);
                    dFBdX[3 * i - 1][3 * i - 4] = 1 / Math.Pow(2, 0.5) * pi * 0.25 * P / nr * ((x[i - 1] - x[i - 2]) * z[i - 2] / r[i - 2]);
                    dFBdX[3 * i - 1][3 * i - 3] = 1 / Math.Pow(2, 0.5) * pi * 0.25 * P / nr * (r[i - 1] + r[i - 2]);
                    dFBdX[3 * i - 1][3 * i - 2] = 0;
                    dFBdX[3 * i - 1][3 * i - 1] = 0;
                }
                else                   // pressure backward starting from npp+1
                {
                    // Fx with respect to x y z (this and previous node)
                    dFBdX[3 * i - 3][3 * i - 6] = 0;
                    dFBdX[3 * i - 3][3 * i - 5] = pi * 0.25 * P / nr * 2 * y[i - 2];
                    dFBdX[3 * i - 3][3 * i - 4] = pi * 0.25 * P / nr * 2 * z[i - 2];
                    dFBdX[3 * i - 3][3 * i - 3] = 0;
                    dFBdX[3 * i - 3][3 * i - 2] = pi * 0.25 * P / nr * (-2 * y[i - 1]);
                    dFBdX[3 * i - 3][3 * i - 1] = pi * 0.25 * P / nr * (-2 * z[i - 1]);
                    // Fy with respect to x y z (this and previous node)
                    dFBdX[3 * i - 2][3 * i - 6] = pi * 0.25 * P / nr * (-y[i - 1] * (r[i - 2] + r[i - 1]) / r[i - 1]);
                    dFBdX[3 * i - 2][3 * i - 5] = pi * 0.25 * P / nr * ((x[i - 1] - x[i - 2]) * y[i - 1] * y[i - 2] / r[i - 1] / r[i - 2]);
                    dFBdX[3 * i - 2][3 * i - 4] = pi * 0.25 * P / nr * ((x[i - 1] - x[i - 2]) * y[i - 1] * z[i - 2] / r[i - 1] / r[i - 2]);
                    dFBdX[3 * i - 2][3 * i - 3] = pi * 0.25 * P / nr * (y[i - 1] * (r[i - 2] + r[i - 1]) / r[i - 1]);
                    dFBdX[3 * i - 2][3 * i - 2] = pi * 0.25 * P / nr * ((x[i - 1] - x[i - 2]) * (r[i - 2] + r[i - 1]) / r[i - 1] - (x[i - 1] - x[i - 2]) * Math.Pow(y[i - 1], 2) * (r[i - 2] + r[i - 1]) / Math.Pow(r[i - 1], 3) + (x[i - 1] - x[i - 2]) * Math.Pow(y[i - 1], 2) / Math.Pow(r[i - 1], 2));
                    dFBdX[3 * i - 2][3 * i - 1] = pi * 0.25 * P / nr * ((x[i - 1] - x[i - 2]) * y[i - 1] * z[i - 1] / Math.Pow(r[i - 1], 2) - (x[i - 1] - x[i - 2]) * y[i - 1] * z[i - 1] * (r[i - 2] + r[i - 1]) / Math.Pow(r[i - 1], 3));
                    // Fz with respect to x y z (this and previous node)
                    dFBdX[3 * i - 1][3 * i - 6] = pi * 0.25 * P / nr * (-z[i - 1] * (r[i - 2] + r[i - 1]) / r[i - 1]);
                    dFBdX[3 * i - 1][3 * i - 5] = pi * 0.25 * P / nr * ((x[i - 1] - x[i - 2]) * z[i - 1] * y[i - 2] / r[i - 1] / r[i - 2]);
                    dFBdX[3 * i - 1][3 * i - 4] = pi * 0.25 * P / nr * ((x[i - 1] - x[i - 2]) * z[i - 1] * z[i - 2] / r[i - 1] / r[i - 2]);
                    dFBdX[3 * i - 1][3 * i - 3] = pi * 0.25 * P / nr * (z[i - 1] * (r[i - 2] + r[i - 1]) / r[i - 1]);
                    dFBdX[3 * i - 1][3 * i - 2] = pi * 0.25 * P / nr * ((x[i - 1] - x[i - 2]) * y[i - 1] * z[i - 1] / Math.Pow(r[i - 1], 2) - (x[i - 1] - x[i - 2]) * y[i - 1] * z[i - 1] * (r[i - 2] + r[i - 1]) / Math.Pow(r[i - 1], 3));
                    dFBdX[3 * i - 1][3 * i - 1] = pi * 0.25 * P / nr * ((x[i - 1] - x[i - 2]) * (r[i - 2] + r[i - 1]) / r[i - 1] - (x[i - 1] - x[i - 2]) * Math.Pow(z[i - 1], 2) * (r[i - 2] + r[i - 1]) / Math.Pow(r[i - 1], 3) + (x[i - 1] - x[i - 2]) * Math.Pow(z[i - 1], 2) / Math.Pow(r[i - 1], 2));
                }
            }
            #endregion

            //TENSION DUE TO CATCH PRESSURE FORWARD
            #region
            for (int i = ncp; i <= 2 * nx; i++)
            {
                if (i == 2 * nx)   // pre last point
                {
                    // Fx with respect to x y z (this and next node)
                    dFFdX[3 * i - 3][3 * i - 3] = pi * 0.25 * P / nr * (0);
                    dFFdX[3 * i - 3][3 * i - 2] = pi * 0.25 * P / nr * (2 * y[i - 1]);
                    dFFdX[3 * i - 3][3 * i - 1] = pi * 0.25 * P / nr * (2 * z[i - 1]);
                    dFFdX[3 * i - 3][3 * i - 0] = pi * 0.25 * P / nr * (0);
                    dFFdX[3 * i - 3][3 * i + 1] = pi * 0.25 * P / nr * (-2 * y[i]);
                    dFFdX[3 * i - 3][3 * i + 2] = pi * 0.25 * P / nr * (-2 * z[i]);
                    // Fy with respect to x y z (this and next node)
                    dFFdX[3 * i - 2][3 * i - 3] = pi * 0.25 * P / nr * (-y[i - 1] * (r[i - 1] + r[i]) / r[i - 1]);
                    dFFdX[3 * i - 2][3 * i - 2] = pi * 0.25 * P / nr * ((x[i] - x[i - 1]) * (r[i] + r[i - 1]) / r[i - 1] - (x[i] - x[i - 1]) * Math.Pow(y[i - 1], 2) * (r[i] + r[i - 1]) / Math.Pow(r[i - 1], 3) + (x[i] - x[i - 1]) * Math.Pow(y[i - 1], 2) / Math.Pow(r[i - 1], 2));
                    dFFdX[3 * i - 2][3 * i - 1] = pi * 0.25 * P / nr * ((x[i] - x[i - 1]) * y[i - 1] * z[i - 1] / Math.Pow(r[i - 1], 2) - (x[i] - x[i - 1]) * y[i - 1] * z[i - 1] * (r[i] + r[i - 1]) / Math.Pow(r[i - 1], 3));
                    dFFdX[3 * i - 2][3 * i - 0] = pi * 0.25 * P / nr * (y[i - 1] * (r[i - 1] + r[i]) / r[i - 1]);
                    dFFdX[3 * i - 2][3 * i + 1] = pi * 0.25 * P / nr * ((x[i] - x[i - 1]) * y[i - 1] / r[i - 1]);
                    dFFdX[3 * i - 2][3 * i + 2] = pi * 0.25 * P / nr * ((x[i] - x[i - 1]) * y[i - 1] / r[i - 1]);
                    // Fz with respect to x y z (this and next node)
                    dFFdX[3 * i - 1][3 * i - 3] = pi * 0.25 * P / nr * (-z[i - 1] * (r[i - 1] + r[i]) / r[i - 1]);
                    dFFdX[3 * i - 1][3 * i - 2] = pi * 0.25 * P / nr * ((x[i] - x[i - 1]) * y[i - 1] * z[i - 1] / Math.Pow(r[i - 1], 2) - (x[i] - x[i - 1]) * y[i - 1] * z[i - 1] * (r[i] + r[i - 1]) / Math.Pow(r[i - 1], 3));
                    dFFdX[3 * i - 1][3 * i - 1] = pi * 0.25 * P / nr * ((x[i] - x[i - 1]) * (r[i] + r[i - 1]) / r[i - 1] - (x[i] - x[i - 1]) * Math.Pow(z[i - 1], 2) * (r[i] + r[i - 1]) / Math.Pow(r[i - 1], 3) + (x[i] - x[i - 1]) * Math.Pow(z[i - 1], 2) / Math.Pow(r[i - 1], 2));
                    dFFdX[3 * i - 1][3 * i - 0] = pi * 0.25 * P / nr * (z[i - 1] * (r[i - 1] + r[i]) / r[i - 1]);
                    dFFdX[3 * i - 1][3 * i + 1] = pi * 0.25 * P / nr * ((x[i] - x[i - 1]) * z[i - 1] / r[i - 1]);
                    dFFdX[3 * i - 1][3 * i + 2] = pi * 0.25 * P / nr * ((x[i] - x[i - 1]) * z[i - 1] / r[i - 1]);
                }
                else        // pressure backward starting from npp+1
                {
                    // Fx with respect to x y z (this and next node)
                    dFFdX[3 * i - 3][3 * i - 3] = pi * 0.25 * P / nr * (0);
                    dFFdX[3 * i - 3][3 * i - 2] = pi * 0.25 * P / nr * (2 * y[i - 1]);
                    dFFdX[3 * i - 3][3 * i - 1] = pi * 0.25 * P / nr * (2 * z[i - 1]);
                    dFFdX[3 * i - 3][3 * i - 0] = pi * 0.25 * P / nr * (0);
                    dFFdX[3 * i - 3][3 * i + 1] = pi * 0.25 * P / nr * (-2 * y[i]);
                    dFFdX[3 * i - 3][3 * i + 2] = pi * 0.25 * P / nr * (-2 * z[i]);
                    // Fy with respect to x y z (this and next node)
                    dFFdX[3 * i - 2][3 * i - 3] = pi * 0.25 * P / nr * (-y[i - 1] * (r[i - 1] + r[i]) / r[i - 1]);
                    dFFdX[3 * i - 2][3 * i - 2] = pi * 0.25 * P / nr * ((x[i] - x[i - 1]) * (r[i] + r[i - 1]) / r[i - 1] - (x[i] - x[i - 1]) * Math.Pow(y[i - 1], 2) * (r[i] + r[i - 1]) / Math.Pow((r[i - 1]), 3) + (x[i] - x[i - 1]) * Math.Pow(y[i - 1], 2) / Math.Pow((r[i - 1]), 2));
                    dFFdX[3 * i - 2][3 * i - 1] = pi * 0.25 * P / nr * ((x[i] - x[i - 1]) * y[i - 1] * z[i - 1] / Math.Pow((r[i - 1]), 2) - (x[i] - x[i - 1]) * y[i - 1] * z[i - 1] * (r[i] + r[i - 1]) / Math.Pow((r[i - 1]), 3));
                    dFFdX[3 * i - 2][3 * i - 0] = pi * 0.25 * P / nr * (y[i - 1] * (r[i - 1] + r[i]) / r[i - 1]);
                    dFFdX[3 * i - 2][3 * i + 1] = pi * 0.25 * P / nr * ((x[i] - x[i - 1]) * y[i] * y[i - 1] / r[i] / r[i - 1]);
                    dFFdX[3 * i - 2][3 * i + 2] = pi * 0.25 * P / nr * ((x[i] - x[i - 1]) * y[i - 1] * z[i] / r[i] / r[i - 1]);
                    // Fz with respect to x y z (this and next node)
                    dFFdX[3 * i - 1][3 * i - 3] = pi * 0.25 * P / nr * (-z[i - 1] * (r[i - 1] + r[i]) / r[i - 1]);
                    dFFdX[3 * i - 1][3 * i - 2] = pi * 0.25 * P / nr * ((x[i] - x[i - 1]) * y[i - 1] * z[i - 1] / Math.Pow((r[i - 1]), 2) - (x[i] - x[i - 1]) * y[i - 1] * z[i - 1] * (r[i] + r[i - 1]) / Math.Pow((r[i - 1]), 3));
                    dFFdX[3 * i - 1][3 * i - 1] = pi * 0.25 * P / nr * ((x[i] - x[i - 1]) * (r[i] + r[i - 1]) / r[i - 1] - (x[i] - x[i - 1]) * Math.Pow(z[i - 1], 2) * (r[i] + r[i - 1]) / Math.Pow((r[i - 1]), 3) + (x[i] - x[i - 1]) * Math.Pow(z[i - 1], 2) / Math.Pow((r[i - 1]), 2));
                    dFFdX[3 * i - 1][3 * i - 0] = pi * 0.25 * P / nr * (z[i - 1] * (r[i - 1] + r[i]) / r[i - 1]);
                    dFFdX[3 * i - 1][3 * i + 1] = pi * 0.25 * P / nr * ((x[i] - x[i - 1]) * y[i] * z[i - 1] / r[i] / r[i - 1]);
                    dFFdX[3 * i - 1][3 * i + 2] = pi * 0.25 * P / nr * ((x[i] - x[i - 1]) * z[i - 1] * z[i] / r[i] / r[i - 1]);
                }
            }
            #endregion

            //dFCdX = dFFdX.Add(dFBdX); // dFCdX = dFFdX + dFBdX collect all contributions at each node
            //Print(dFFdX, new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 }, new int[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 });
            BandAdd2(ref dFCdX, dFFdX, dFBdX);
            return dFCdX;
        }

        public double[][] GetTotalJacobian(double[] X)
        {
            dFTdX = GetTwineJacobian(X);
            dFCdX = GetCatchJacobian(X);

            //dFdX = dFTdX.Add(dFCdX);            /*Add elastic twine jacobian and catch presure jacobian*/

            BandAdd2(ref dFdX, dFTdX, dFCdX);

            return dFdX;
        }

        /*check numaerical jacobian*/

        public void CheckJacobian(double[][] J, double[] X)
        {
            double h = Math.Pow(10, -9);
            double[] Xpert = new double[X.GetLength(0)];
            double[] F, Fpert;
            double[][] numJ = new double[X.GetLength(0)][];
            for (int i = 0; i < X.GetLength(0); i++)
            {
                numJ[i] = new double[X.GetLength(0)];
            }
            List<string> Mismatch = new List<string>();
            for (int i = 0; i < X.GetLength(0); i++)
            {
                X.CopyTo(Xpert);
                Xpert[i] += h;

                F = GetTotalForces(X);
                Fpert = GetTotalForces(Xpert);

                for (int j = 0; j < X.GetLength(0); j++)
                {
                    numJ[j][i] = (Fpert[j] - F[j]) / h - J[j][i];

                    if (numJ[j][i] > Math.Pow(10, -2))
                    {
                        Mismatch.Add(string.Format("{0};{1}", i, j));
                    }
                }
            }
        }

            /* printing*/

        public void Print(double[] Array)
        {
            Console.WriteLine();

            for (int i = 0; i < Array.Length; i++)
            {
                Console.WriteLine(String.Format("{0}", Array[i]));
            }
        }

        public void Print(double[][] Array, int[] Rows, int[] Cols)
        {
            Console.WriteLine();
            for (int i = 0; i < Rows.Length; i++)
            {
                for (int j = 0; j < Cols.Length; j++)
                {
                    Console.Write(String.Format("{0,-15:E4}", Array[Rows[i]][Cols[j]]));
                }
                Console.WriteLine();
            }
        }

        public void PrintResults(double[] X)
        {
            for (int i = 0; i < dof / 3; i++)
            {
                r[i] = Math.Sqrt(Math.Pow(X[3 * i + 1], 2) + Math.Pow(X[3 * i + 2], 2));
            }

            /*Maximum radius and total length*/
            #region
            double L = X[dof - 3];                      // total length is last x
            double rCatch = r[ncp - 1];                 // radius at the beginning of the catch
            double rMax = r0;                           // maximum radius

            for (int i = 1; i < dof / 3; i++)
            {
                if (rMax < r[i])
                {
                    rMax = r[i];
                }
            }
            #endregion

            /* Get beginnig reaction and drag force on the catch*/
            #region

            F = GetTotalForces(X);
            double begR = 2 * nr * F[0];                                   // resultant entrance reaction
            double dragR = pi * Math.Pow(rCatch, 2) * P;               // resultant catch drag

            #endregion

            /*Length, surface area and volume covered by the catch*/
            #region
            double Lc = X[dof - 3] - X[3 * ncp - 3];    // catch length difference between last x and x where catch starts

            double S = 0;                               // surface area affected by catch
            double V = 0;                               // volume of the catch

            for (int i = ncp; i < dof / 3; i++)         // surface area of a cut cone between 2 nodes
            {                                           // volume of a cut cone between 2 nodes
                S = S + pi * (r[i] + r[i-1]) * Math.Pow(Math.Pow(r[i-1] - r[i], 2) + Math.Pow(X[3 * i] - X[3 * i - 3], 2), 0.5);
                V = V + pi * (X[3 * i] - X[3 * i - 3]) * (Math.Pow(r[i], 2) + r[i] * r[i - 1] + Math.Pow(r[i - 1], 2)) / 3;
            }

            #endregion

            /*Print all in a table*/
            #region
            Console.WriteLine("\nResults:");
            Console.WriteLine("Total codend length:             {0,10:N3} [m]", L);
            Console.WriteLine("Maximum codend radius:           {0,10:N3} [m]", rMax);
            Console.WriteLine("Catch extent:                    {0,10:N3} [m]", Lc);
            Console.WriteLine("Surface in contact with catch:   {0,10:N3} [m^2]", S);
            Console.WriteLine("Catch volume:                    {0,10:N3} [m^3]", V);
            Console.WriteLine("Resultant entrance reaction:     {0,10:N0} [N]", begR);
            Console.WriteLine("Total drag force:                {0,10:N0} [N]\n", dragR);
            #endregion

        }

        public void PrintIter(int iter, int funEval, double addStiff, double R)
        {
            Console.WriteLine(String.Format("Iteration: {0,-10:D} " +
                                "Line search steps: {1,-10:D} " +
                                "Added stiffness: {2,-10:N3} " +
                                "Force residue: {3,-10:e3}", iter, funEval, addStiff, R));
        }

    }
}
