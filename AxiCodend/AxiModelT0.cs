using System;
using System.Collections.Generic;
using System.IO;
using CSparse.Storage;
using CSparse;

namespace AxiCodend
{
    class AxiModelT0
    {
        // CLASS VARIABLES
        #region
        /*codend paramters*/
        public int nx;              // amount of meshes along the codend length
        public int nr;              // amount of meshes in the codend`s circumference
        public double r0;           // entrance radius
        public int nc;              // moundt of meshes along the codend length contatining catch 

        /*mesh characteristics*/
        private double l0;           // length of un-stretched knot       
        private double m0;           // length of un-stretched twine
        private double kl;           // stiffness of a knot
        private double km;           // stiffness of a twine

        /* operation */
        public double towSpeed;
        private double P;            // catch pressure
        #endregion

        public HexMeshPanelMaterial Material;
        /* calculation variables*/
        #region
        public readonly int dof;            // nb of dof
        public readonly int nbc_beg = 2;	// nb of fixed dof at entry (x0 r0)
        public readonly int nbc_end = 1;    // nb of fixed dof at the end (rn)

        const double pi = Math.PI;
        private int ncp;                          // node number indicating start of application of the catch pressure
        private double theta, co, si, co2, si2;   // angle between 2 radial planes + sin cos of single and double theta

        public double R;                         // residue
        public double[] X;                       // Shape
        public double[] F;                       // force vectors
        private CoordinateStorage<double> Jcs;   // Sparse Jacobians
        public CompressedColumnStorage<double> J;

        /*help figures and help variables*/
        #region
        double a1, a2, a3, a4, a5, a6, l10, l11;

        // neigbour points (local dof) to the left and to the right of "a"
        //					
        //          (a5,a6) _______
        //                 /       \
        //                /l11      \
        //       ________/           \________     -> current
        // (a1,a2) l10  (a3,a4) 
        //
        double b1, b2, b3, b4, b5, b6, l21, l22;

        // neigbour points (local dof) to the left and to the right of "b"
        //					
        //           (b3,b4)_______(b5,b6) 
        //                 /  l22  \
        //                /l21      \
        //       ________/           \________     -> current
        //			  (b1,b2) 
        //
        double c1, c2, c3, c4, c5, c6, l33;

        // neigbour points (local dof) to the left and to the right of "c"
        //					
        //			 (c1,c2)_______(c3,c4) 
        //                 /	   \
        //                /		 l33\
        //       ________/           \________     -> current
        //						(c5,c6) 
        //
        double d1, d2, d3, d4, d5, d6, l44;

        // neigbour points (local dof) to the left and to the right of "d"
        //					
        //					_______(d1,d2) 
        //                 /	   \
        //                /			\
        //       ________/           \________     -> current
        //						(d3,d4)	 l44  (d5,d6) 
        //
        #endregion

        #endregion

        // CLASS CONSTRUCTS

        public AxiModelT0(int nx, int nr, double r0, HexMeshPanelMaterial Material)
        {
            /*main inputs*/
            this.nx = nx;
            this.nr = nr;
            this.r0 = r0;
            this.Material = Material;

            l0 = Material.KnotSize;
            m0 = Material.MeshSide / 2;
            kl = Material.KnotEA;
            km = Material.TwineEA;

            dof = nx * 8 + 2;

            // default towing and catch to start with
            nc = nx;
            ncp = 1 + 4 * nx - 4 * nc;

            towSpeed = 1;
            P = 0.5 * 1.4 * 1025.0 * Math.Pow(towSpeed, 2);

            /*initialize state*/
            ClearState();

            /*angles and trigonometric functions*/
            SetAngles();
        }

        public AxiModelT0(PathsIO path)
        {
            /*main inputs*/
            LoadInput(path);

            dof = nx * 8 + 2;

            // default towing and catch to start with
            nc = nx;
            ncp = 1 + 4 * nx - 4 * nc;

            towSpeed = 1;
            P = 0.5 * 1.4 * 1025.0 * Math.Pow(towSpeed, 2);

            /*initialize state*/
            ClearState();

            /*angles and trigonometric functions*/
            SetAngles();
        }


        // CLASS METHODS

        private void LoadInput(PathsIO path)
        {
            string[] names = { "MeshSide", "KnotSize", "TwineEA", "KnotEA",
                               "MeshesAlong", "MeshesAround", "EntranceRadius" };
            string[] lines = File.ReadAllLines(path.input);
            string[] parts;
            int currentLine = 0;
            int foundCount = 0;

            Material = new HexMeshPanelMaterial(path);

            foreach (var line in lines)
            {
                if (line.Contains(names[0]))
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    m0 = Convert.ToDouble(parts[1]) / 2;
                    foundCount++;
                }

                if (line.Contains(names[1]))
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    l0 = Convert.ToDouble(parts[1]);
                    foundCount++;
                }

                if (line.Contains(names[2]))
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    km = Convert.ToDouble(parts[1]);
                    foundCount++;
                }

                if (line.Contains(names[3]))
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    kl = Convert.ToDouble(parts[1]);
                    foundCount++;
                }

                if (line.Contains(names[4]))
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    nx = Convert.ToInt32(parts[1]);
                    foundCount++;
                }

                if (line.Contains(names[5]))
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    nr = Convert.ToInt32(parts[1]);
                    foundCount++;
                }

                if (line.Contains(names[6]))
                {
                    parts = lines[currentLine].Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    r0 = Convert.ToDouble(parts[1]);
                    foundCount++;
                }
                currentLine++;
            }

            if (foundCount != 7)
            {
                throw new ArgumentException("Not all fields could be initialized, " +
                                            "because the input file is not in the right format");
            }
        }

        /* for updating*/

        private void SetAngles()
        {
            theta = pi / nr;
            co = Math.Cos(theta);
            si = Math.Sin(theta);
            co2 = Math.Cos(2 * theta);
            si2 = Math.Sin(2 * theta);
        }

        public void ClearState()
        {
            X = new double[dof];
            F = new double[dof];
            R = 0;
            Jcs = new CoordinateStorage<double>(dof, dof, 4 * 4 + (dof - 4) * 6);
        }

        public void UpdatePosition(double[] h, double lambda)
        {
            for (int i = nbc_beg; i < dof - nbc_end; i++)
            {
                X[i] = X[i] - lambda * h[i];
            }
        }

        public void ApplyCatch(int newCatch)
        {
            nc = newCatch;
            ncp = 1 + 4 * nx - 4 * nc;
        }

        public void ApplyTowing(double newSpeed)
        {
            towSpeed = newSpeed;
            P = 0.5 * 1.4 * 1025.0 * Math.Pow(towSpeed, 2);
        }

        /* initial shape*/

        public void SetInitialShape()
        {
            for (int i = 0; i < dof / 2; i++)
            {
                X[2 * i] = 0.5 * (l0 + m0) * i;
                X[2 * i + 1] = r0;
            }

            X[dof - 1] = 0;
        }

        public void SetInitialShapeSmooth()
        {
            double fStretch = 1.005;
            int Lint = (int)(pi * r0 / (fStretch * (l0 + m0)));

            for (int i = 0; i < dof / 2 - Lint; i++)
            {
                X[2 * i + 1] = r0;

                if (i == 0)
                {
                    continue;
                }

                if (i % 2 == 1)
                {
                    X[2 * i] = X[2 * i - 2] + fStretch * l0;
                }
                else
                {
                    X[2 * i] = X[2 * i - 2] + fStretch * m0;
                }
                //Console.WriteLine(i + "\t" + X0[2 * i] + "\t" +X0[2 * i + 1]);
            }

            double omega = 0.5 * pi;

            //Console.WriteLine(dof - 2 * L - 2 + " " + X0[dof - 2 * Lint - 2]);
            for (int i = 0; i < Lint; i++)
            {
                if (i % 2 == 0)
                {
                    omega -= 0.5 * pi * fStretch * l0 / (pi * r0 / 2);
                }
                else
                {
                    omega -= 0.5 * pi * fStretch * m0 / (pi * r0 / 2);
                }

                X[2 * i + dof - 2 * Lint] = X[dof - 2 * Lint - 2] + r0 * Math.Cos(omega);
                X[2 * i + 1 + dof - 2 * Lint] = r0 * Math.Sin(omega);
                //Console.WriteLine(2 * i + dof - 2 * Lint + "\t" + omega +"\t" + X0[2 * i + dof - 2 * Lint] + "\t" + X0[2 * i + 1 + dof - 2 * Lint]);
            }
            X[dof - 1] = 0;
        }

        public void SetInitialShape(string path)
        {
            string[] lines = File.ReadAllLines(path);

            for (int i = 0; i < lines.Length; i++)
            {
                X[i] = Convert.ToDouble(lines[i]);
            }
        }

        /* forces*/

        private void SetTwineForces(bool IncludeBC)
        {
            /*Calculate forces at each node a-b-c-d*/

            for (int i = 0; i < nx; i++)
            {
                /* node a */
                #region
                a1 = X[8 * i + 0];  // dx (i-1)
                a2 = X[8 * i + 1];  // dr (i-1)
                a3 = X[8 * i + 2];  // ax (i)
                a4 = X[8 * i + 3];  // ar (i)
                a5 = X[8 * i + 4];  // bx (i)
                a6 = X[8 * i + 5];  // br (i)

                l10 = Math.Sqrt(Math.Pow((a3 - a1), 2) + Math.Pow((a4 - a2), 2));                      // d-a knot length
                if (l10 < l0)
                    l10 = l0;

                F[8 * i + 2] += kl * ((l10 - l0)) / l0 / l10 * (a1 - a3);                             // x component due to tension in d-a knot
                F[8 * i + 3] += kl * ((l10 - l0)) / l0 / l10 * (a2 - a4);                             // r component due to tension in d-a knot

                l11 = Math.Sqrt(Math.Pow((a5 - a3), 2) + Math.Pow((a6 * si), 2) + Math.Pow((a6 * co - a4), 2));             // a-b twine length
                if (l11 < m0)
                    l11 = m0;

                F[8 * i + 2] += 2 * km * (l11 - m0) / m0 / l11 * (a5 - a3);            // x component due to tension in a-b twine (above and below)
                F[8 * i + 3] += 2 * km * (l11 - m0) / m0 / l11 * (a6 * co - a4);       // r component due to tension in a-b twine (above and below)
                #endregion

                /* node b */
                #region
                b1 = X[8 * i + 2]; // ax (i)
                b2 = X[8 * i + 3]; // ar (i)
                b3 = X[8 * i + 4]; // bx (i)
                b4 = X[8 * i + 5]; // br (i)
                b5 = X[8 * i + 6]; // cx (i)
                b6 = X[8 * i + 7]; // cr (i)

                l21 = Math.Sqrt(Math.Pow((b3 - b1), 2) + Math.Pow((b4 * si), 2) + Math.Pow((b4 * co - b2), 2));                  // a-b twine length
                if (l21 < m0)
                    l21 = m0;

                F[8 * i + 4] += 2 * km * (l21 - m0) / m0 / l21 * (b1 - b3);                                                      // x component due to tension in a-b twine (above and below)                                                                                 
                F[8 * i + 5] += 2 * km * (l21 - m0) / m0 / l21 * (b2 * co - b4);                                                 // r component due to tension in a-b twine (above and below)

                l22 = Math.Sqrt(Math.Pow((b5 - b3), 2) + Math.Pow((b6 * si - b4 * si), 2) + Math.Pow((b6 * co - b4 * co), 2));   // b-c knot length
                if (l22 < l0)
                    l22 = l0;

                F[8 * i + 4] += kl * (l22 - l0) / l0 / l22 * (b5 - b3);                                          // x component due to tension in b-c knot
                F[8 * i + 5] += kl * (l22 - l0) / l0 / l22 * (b6 - b4);                                          // r component due to tension in b-c knot
                #endregion

                /* node c */
                #region
                c1 = X[8 * i + 4]; // bx (i)
                c2 = X[8 * i + 5]; // br (i)
                c3 = X[8 * i + 6]; // cx (i)
                c4 = X[8 * i + 7]; // cr (i)
                c5 = X[8 * i + 8]; // dx (i)
                c6 = X[8 * i + 9]; // dr (i)

                l22 = Math.Sqrt(Math.Pow((c3 - c1), 2) + Math.Pow((c4 * si - c2 * si), 2) + Math.Pow((c4 * co - c2 * co), 2));  // b-c knot length
                if (l22 < l0)
                    l22 = l0;

                F[8 * i + 6] += kl * (l22 - l0) / l0 / l22 * (c1 - c3);                                         // x component due to tension in b-c knot
                F[8 * i + 7] += kl * (l22 - l0) / l0 / l22 * (c2 - c4);                                         // r component due to tension in b-c knot

                l33 = Math.Sqrt(Math.Pow((c5 - c3), 2) + Math.Pow((-c4 * si), 2) + Math.Pow((c6 - c4 * co), 2));    // c-d twine length
                if (l33 < m0)
                    l33 = m0;

                F[8 * i + 6] += 2 * km * (l33 - m0) / l33 / m0 * (c5 - c3);                    // x component due to tension in c-d knot
                F[8 * i + 7] += 2 * km * (l33 - m0) / l33 / m0 * (c6 * co - c4);               // r component due to tension in c-d knot

                #endregion

                /*node d (first mesh)*/
                #region
                if (i == 0)
                {
                    if (IncludeBC)
                    {
                        F[0] = F[1] = 0;
                    }
                    else
                    {
                        d3 = X[0]; // dx(0)
                        d4 = X[1]; // dr(0)
                        d5 = X[2]; // ax(1)
                        d6 = X[3]; // ar(1)

                        l44 = Math.Sqrt(Math.Pow((d5 - d3), 2) + Math.Pow((d6 - d4), 2));          // d-a knot length
                        if (l44 < l0)
                            l44 = l0;

                        F[0] += kl * (l44 - l0) / l44 / l0 * (d5 - d3);                             // x component due to tension in d-a knot
                        F[1] += kl * (l44 - l0) / l44 / l0 * (d6 - d4);                            // r component due to tension in d-a knot
                    }
                }
                #endregion

                /*node d (last mesh)*/
                #region
                if (i == nx - 1)
                {
                    d1 = X[8 * i + 6]; // cx(last)
                    d2 = X[8 * i + 7]; // cr(last)
                    d3 = X[8 * i + 8]; // dx(last)
                    d4 = X[8 * i + 9]; // dr(last)

                    l33 = Math.Sqrt(Math.Pow((d3 - d1), 2) + Math.Pow((-d2 * si), 2) + Math.Pow((d4 - d2 * co), 2));                // c-d knot length
                    if (l33 < m0)
                        l33 = m0;

                    F[8 * i + 8] += -km * (l33 - m0) / l33 / m0 * (d3 - d1);                                   // x component due to tension in c-d knot (above and below)
                    F[8 * i + 9] += -km * (l33 - m0) / l33 / m0 * (d4 - d2 * co);                              // x component due to tension in c-d knot (above and below)

                    F[8 * i + 8] += km * (l33 - m0) / l33 / m0 * (d1 - d3);
                    F[8 * i + 9] += km * (l33 - m0) / l33 / m0 * (d2 * co - d4);

                    if (IncludeBC)
                    {
                        F[8 * i + 9] = 0;
                    }
                }
                #endregion

                /* node d (intermediate meshes)*/
                #region
                else
                {
                    d1 = X[8 * i + 6]; // cx (i)
                    d2 = X[8 * i + 7]; // cr (i)
                    d3 = X[8 * i + 8]; // dx (i)
                    d4 = X[8 * i + 9]; // dr (i)
                    d5 = X[8 * i + 10]; // ax (i+1)
                    d6 = X[8 * i + 11]; // ar (i+1)

                    l33 = Math.Sqrt(Math.Pow((d3 - d1), 2) + Math.Pow((-d2 * si), 2) + Math.Pow((d4 - d2 * co), 2));        // c-d twine length
                    if (l33 < m0)
                        l33 = m0;

                    F[8 * i + 8] += 2 * km * (l33 - m0) / l33 / m0 * (d1 - d3);                                          // x component due to tension in c-d twine
                    F[8 * i + 9] += 2 * km * (l33 - m0) / l33 / m0 * (d2 * co - d4);                                     // r component due to tension in c-d twine

                    l44 = Math.Sqrt(Math.Pow((d5 - d3), 2) + Math.Pow((d6 - d4), 2));                                   // d-a knot length
                    if (l44 < l0)
                        l44 = l0;

                    F[8 * i + 8] += kl * (l44 - l0) / l44 / l0 * (d5 - d3);                            // x component due to tension in d-a knot
                    F[8 * i + 9] += kl * (l44 - l0) / l44 / l0 * (d6 - d4);                            // r component due to tension in d-a knot
                    #endregion
                }
            }
        }

        private void SetCatchForces(bool IncludeBC)
        {
            /* intermediate meshes */
            #region
            for (int i = ncp + 1; i < dof / 2; i++)
            {
                F[2 * i - 2] += 0.5 * P * theta * (Math.Pow(X[2 * i - 3], 2) - Math.Pow(X[2 * i + 1], 2));                  // Axial component
                F[2 * i - 1] += 0.5 * P * theta * (X[2 * i - 1] + X[2 * i + 1]) * Math.Abs(X[2 * i] - X[2 * i - 2]) +       // Radial component
                                0.5 * P * theta * (X[2 * i - 3] + X[2 * i - 1]) * Math.Abs(X[2 * i - 2] - X[2 * i - 4]);
            }
            #endregion

            /*first mesh*/
            #region
            F[2 * ncp - 2] += 0.5 * P * theta * (Math.Pow(X[2 * ncp - 1], 2) - Math.Pow(X[2 * ncp + 1], 2));                // Axial component
            F[2 * ncp - 1] += 0.5 * P * theta * (X[2 * ncp - 1] + X[2 * ncp + 1]) * Math.Abs(X[2 * ncp] - X[2 * ncp - 2]);  // Radial component

            if (ncp == 1 && IncludeBC)
            {
                F[2 * ncp - 2] = F[2 * ncp - 1] = 0;
            }

            #endregion

            /*last mesh*/
            #region
            F[dof - 2] += 0.5 * P * theta * (Math.Pow(X[dof - 3], 2) - Math.Pow(X[dof - 1], 2));                            // Axial component
            F[dof - 1] += 0.5 * P * theta * (X[dof - 3] + X[dof - 1]) * Math.Abs(X[dof - 2] - X[dof - 4]);                 // Radial component

            if (IncludeBC)
            {
                F[dof - 1] = 0;
            }
            #endregion

        }

        public void UpdateTotalForces(bool IncludeBC)
        {
            Array.Clear(F, 0, F.Length);
            SetTwineForces(IncludeBC);
            SetCatchForces(IncludeBC);
        }

        public void UpdateResidual()
        {
            R = 0;
            for (int i = nbc_beg; i < dof - nbc_end; i++)
            {
                R += Math.Pow(F[i], 2);                      // residue as sum of squares
            }
            R = Math.Sqrt(R);
        }

        /* jacobian / stiffness*/

        private void SetTwineJacobian(bool IncludeBC)
        {
            /*Calculate forces at each node a-b-c-d*/

            for (int i = 0; i < nx; i++)
            {
                /* node a */
                #region
                a1 = X[8 * i + 0];  // dx (i-1)
                a2 = X[8 * i + 1];  // dr (i-1)
                a3 = X[8 * i + 2];  // ax (i)
                a4 = X[8 * i + 3];  // ar (i)
                a5 = X[8 * i + 4];  // bx (i)
                a6 = X[8 * i + 5];  // br (i)

                l10 = Math.Sqrt(Math.Pow((a3 - a1), 2) + Math.Pow((a4 - a2), 2));                                           // d-a knot length
                if (l10 < l0)
                    l10 = l0;

                l11 = Math.Sqrt(Math.Pow((a5 - a3), 2) + Math.Pow((a6 * si), 2) + Math.Pow((a6 * co - a4), 2));             // a-b twine length
                if (l11 < m0)
                    l11 = m0;

                // Stifness matrix elements corresponding to point a - axial direction
                if (i == 0 && IncludeBC)
                { }
                else
                {
                    Jcs.At(8 * i + 2, 8 * i + 0, (kl * (Math.Sqrt(Math.Pow((a4 - a2), 2) + Math.Pow((a3 - a1), 2)) - l0)) / (Math.Sqrt(Math.Pow((a4 - a2), 2) + Math.Pow((a3 - a1), 2)) * l0) - (Math.Pow((a3 - a1), 2) * kl * (Math.Sqrt(Math.Pow((a4 - a2), 2) + Math.Pow((a3 - a1), 2)) - l0)) / (Math.Pow(l10, 3) * l0) + (Math.Pow((a3 - a1), 2) * kl) / (Math.Pow(l10, 2) * l0));
                    Jcs.At(8 * i + 2, 8 * i + 1, (kl * (a4 - a2) * (a3 - a1)) / (l0 * (Math.Pow((a3 - a1), 2) + Math.Pow((a4 - a2), 2))) - (kl * (a4 - a2) * (l10 - l0) * (a3 - a1)) / (l0 * ((Math.Pow(l10, 3)))));
                }
                Jcs.At(8 * i + 2, 8 * i + 2, -(kl * Math.Pow((a3 - a1), 2)) / (l0 * (Math.Pow((a3 - a1), 2) + Math.Pow((a4 - a2), 2))) + (kl * (l10 - l0) * Math.Pow((a3 - a1), 2)) / (l0 * ((Math.Pow(l10, 3)))) - (kl * (l10 - l0)) / (l0 * l10) - (2 * km * Math.Pow((a5 - a3), 2)) / (m0 * Math.Pow(l11, 2)) + (2 * km * (l11 - m0) * Math.Pow((a5 - a3), 2)) / (m0 * Math.Pow(l11, 3)) - (2 * km * (l11 - m0)) / (m0 * l11));
                Jcs.At(8 * i + 2, 8 * i + 3, -(kl * (a4 - a2) * (a3 - a1)) / (l0 * (Math.Pow((a3 - a1), 2) + Math.Pow((a4 - a2), 2))) + (kl * (a4 - a2) * (l10 - l0) * (a3 - a1)) / (l0 * ((Math.Pow(l10, 3)))) - (2 * km * (co * a6 - a4) * (a5 - a3)) / (m0 * Math.Pow(l11, 2)) + (2 * km * (co * a6 - a4) * (l11 - m0) * (a5 - a3)) / (m0 * Math.Pow(l11, 3)));
                Jcs.At(8 * i + 2, 8 * i + 4, (2 * km * Math.Pow((a5 - a3), 2)) / (m0 * Math.Pow(l11, 2)) - (2 * km * (l11 - m0) * Math.Pow((a5 - a3), 2)) / (m0 * Math.Pow(l11, 3)) + (2 * km * (l11 - m0)) / (m0 * l11));
                Jcs.At(8 * i + 2, 8 * i + 5, (km * (2 * co * (co * a6 - a4) + 2 * Math.Pow(si, 2) * a6) * (a5 - a3)) / (m0 * Math.Pow(l11, 2)) - (km * (2 * co * (co * a6 - a4) + 2 * Math.Pow(si, 2) * a6) * (l11 - m0) * (a5 - a3)) / (m0 * Math.Pow(l11, 3)));

                // Stifness matrix elements corresponding to point a - radial direction
                if (i == 0 && IncludeBC)
                { }
                else
                {
                    Jcs.At(8 * i + 3, 8 * i + 0, (kl * (a4 - a2) * (a3 - a1)) / (l0 * (Math.Pow((a3 - a1), 2) + Math.Pow((a4 - a2), 2))) - (kl * (a4 - a2) * (l10 - l0) * (a3 - a1)) / (l0 * ((Math.Pow(l10, 3)))));
                    Jcs.At(8 * i + 3, 8 * i + 1, (kl * (l10 - l0)) / (l0 * l10) + (kl * Math.Pow((a4 - a2), 2)) / (l0 * (Math.Pow((a3 - a1), 2) + Math.Pow((a4 - a2), 2))) - (kl * Math.Pow((a4 - a2), 2) * (l10 - l0)) / (l0 * ((Math.Pow(l10, 3)))));
                }
                Jcs.At(8 * i + 3, 8 * i + 2, -(kl * (a4 - a2) * (a3 - a1)) / (l0 * (Math.Pow((a3 - a1), 2) + Math.Pow((a4 - a2), 2))) + (kl * (a4 - a2) * (l10 - l0) * (a3 - a1)) / (l0 * ((Math.Pow(l10, 3)))) - (2 * km * (co * a6 - a4) * (a5 - a3)) / (m0 * Math.Pow(l11, 2)) + (2 * km * (co * a6 - a4) * (l11 - m0) * (a5 - a3)) / (m0 * Math.Pow(l11, 3)));
                Jcs.At(8 * i + 3, 8 * i + 3, -(kl * (l10 - l0)) / (l0 * l10) - (kl * Math.Pow((a4 - a2), 2)) / (l0 * (Math.Pow((a3 - a1), 2) + Math.Pow((a4 - a2), 2))) + (kl * Math.Pow((a4 - a2), 2) * (l10 - l0)) / (l0 * ((Math.Pow(l10, 3)))) - (2 * km * (l11 - m0)) / (m0 * l11) - (2 * km * Math.Pow((co * a6 - a4), 2)) / (m0 * Math.Pow(l11, 2)) + (2 * km * Math.Pow((co * a6 - a4), 2) * (l11 - m0)) / (m0 * Math.Pow(l11, 3)));
                Jcs.At(8 * i + 3, 8 * i + 4, (2 * km * (co * a6 - a4) * (a5 - a3)) / (m0 * Math.Pow(l11, 2)) - (2 * km * (co * a6 - a4) * (l11 - m0) * (a5 - a3)) / (m0 * Math.Pow(l11, 3)));
                Jcs.At(8 * i + 3, 8 * i + 5, (2 * km * co * (l11 - m0)) / (m0 * l11) + (km * (2 * co * (co * a6 - a4) + 2 * Math.Pow(si, 2) * a6) * (co * a6 - a4)) / (m0 * Math.Pow(l11, 2)) - (km * (2 * co * (co * a6 - a4) + 2 * Math.Pow(si, 2) * a6) * (co * a6 - a4) * (l11 - m0)) / (m0 * Math.Pow(l11, 3)));

                #endregion

                /* node b */
                #region
                b1 = X[8 * i + 2]; // ax (i)
                b2 = X[8 * i + 3]; // ar (i)
                b3 = X[8 * i + 4]; // bx (i)
                b4 = X[8 * i + 5]; // br (i)cos(teta)
                b5 = X[8 * i + 6]; // cx (i)
                b6 = X[8 * i + 7]; // cr (i)

                l21 = Math.Sqrt(Math.Pow((b3 - b1), 2) + Math.Pow((b4 * si), 2) + Math.Pow((b4 * co - b2), 2));                  // a-b twine length
                if (l21 < m0)
                    l21 = m0;

                l22 = Math.Sqrt(Math.Pow((b5 - b3), 2) + Math.Pow((b6 * si - b4 * si), 2) + Math.Pow((b6 * co - b4 * co), 2));   // b-c knot length
                if (l22 < l0)
                    l22 = l0;

                // Stifness matrix elements corresponding to point b - axial direction
                Jcs.At(8 * i + 4, 8 * i + 2, -(km * (b3 - b1) * (b1 - b3)) / (m0 * Math.Pow(l21, 2)) + (km * (l21 - m0) * (b3 - b1) * (b1 - b3)) / (m0 * Math.Pow(l21, 3)) + (km * Math.Pow((b3 - b1), 2)) / (m0 * Math.Pow(l21, 2)) - (km * (l21 - m0) * Math.Pow((b3 - b1), 2)) / (m0 * Math.Pow(l21, 3)) + (2 * km * (l21 - m0)) / (m0 * l21));
                Jcs.At(8 * i + 4, 8 * i + 3, -(km * (co * b4 - b2) * (b1 - b3)) / (m0 * Math.Pow(l21, 2)) + (km * (co * b4 - b2) * (l21 - m0) * (b1 - b3)) / (m0 * Math.Pow(l21, 3)) + (km * (co * b4 - b2) * (b3 - b1)) / (m0 * Math.Pow(l21, 2)) - (km * (co * b4 - b2) * (l21 - m0) * (b3 - b1)) / (m0 * Math.Pow(l21, 3)));
                Jcs.At(8 * i + 4, 8 * i + 4, (km * (b3 - b1) * (b1 - b3)) / (m0 * Math.Pow(l21, 2)) - (km * (l21 - m0) * (b3 - b1) * (b1 - b3)) / (m0 * Math.Pow(l21, 3)) - (km * Math.Pow((b3 - b1), 2)) / (m0 * Math.Pow(l21, 2)) + (km * (l21 - m0) * Math.Pow((b3 - b1), 2)) / (m0 * Math.Pow(l21, 3)) - (2 * km * (l21 - m0)) / (m0 * l21) - (kl * Math.Pow((b5 - b3), 2)) / (l0 * Math.Pow(l22, 2)) + (kl * (l22 - l0) * Math.Pow((b5 - b3), 2)) / (l0 * Math.Pow(l22, 3)) - (kl * (l22 - l0)) / (l0 * l22));
                Jcs.At(8 * i + 4, 8 * i + 5, (km * (2 * co * (co * b4 - b2) + 2 * Math.Pow(si, 2) * b4) * (b1 - b3)) / (2 * m0 * Math.Pow(l21, 2)) - (km * (2 * co * (co * b4 - b2) + 2 * Math.Pow(si, 2) * b4) * (l21 - m0) * (b1 - b3)) / (2 * m0 * Math.Pow(l21, 3)) - (km * (2 * co * (co * b4 - b2) + 2 * Math.Pow(si, 2) * b4) * (b3 - b1)) / (2 * m0 * Math.Pow(l21, 2)) + (km * (2 * co * (co * b4 - b2) + 2 * Math.Pow(si, 2) * b4) * (l21 - m0) * (b3 - b1)) / (2 * m0 * Math.Pow(l21, 3)) + (kl * (-2 * si * (si * b6 - si * b4) - 2 * co * (co * b6 - co * b4)) * (b5 - b3)) / (2 * l0 * Math.Pow(l22, 2)) - (kl * (-2 * si * (si * b6 - si * b4) - 2 * co * (co * b6 - co * b4)) * (l22 - l0) * (b5 - b3)) / (2 * l0 * Math.Pow(l22, 3)));
                Jcs.At(8 * i + 4, 8 * i + 6, (kl * Math.Pow((b5 - b3), 2)) / (l0 * Math.Pow(l22, 2)) - (kl * (l22 - l0) * Math.Pow((b5 - b3), 2)) / (l0 * Math.Pow(l22, 3)) + (kl * (l22 - l0)) / (l0 * l22));
                Jcs.At(8 * i + 4, 8 * i + 7, (kl * (2 * si * (si * b6 - si * b4) + 2 * co * (co * b6 - co * b4)) * (b5 - b3)) / (2 * l0 * Math.Pow(l22, 2)) - (kl * (2 * si * (si * b6 - si * b4) + 2 * co * (co * b6 - co * b4)) * (l22 - l0) * (b5 - b3)) / (2 * l0 * Math.Pow(l22, 3)));

                // Stifness matrix elements corresponding to point b - radial direction
                Jcs.At(8 * i + 5, 8 * i + 2, -(km * (co * b2 - b4) * (b3 - b1)) / (m0 * Math.Pow(l21, 2)) + (km * (b4 - co * b2) * (b3 - b1)) / (m0 * Math.Pow(l21, 2)) + (km * (co * b2 - b4) * (l21 - m0) * (b3 - b1)) / (m0 * Math.Pow(l21, 3)) - (km * (b4 - co * b2) * (l21 - m0) * (b3 - b1)) / (m0 * Math.Pow(l21, 3)));
                Jcs.At(8 * i + 5, 8 * i + 3, (km * co * (l21 - m0)) / (m0 * l21) + (km * co * (l21 - m0)) / (m0 * l21) - (km * (co * b4 - b2) * (co * b2 - b4)) / (m0 * Math.Pow(l21, 2)) + (km * (co * b4 - b2) * (b4 - co * b2)) / (m0 * Math.Pow(l21, 2)) + (km * (co * b4 - b2) * (co * b2 - b4) * (l21 - m0)) / (m0 * Math.Pow(l21, 3)) - (km * (co * b4 - b2) * (b4 - co * b2) * (l21 - m0)) / (m0 * Math.Pow(l21, 3)));
                Jcs.At(8 * i + 5, 8 * i + 4, (km * (co * b2 - b4) * (b3 - b1)) / (m0 * Math.Pow(l21, 2)) - (km * (b4 - co * b2) * (b3 - b1)) / (m0 * Math.Pow(l21, 2)) - (km * (co * b2 - b4) * (l21 - m0) * (b3 - b1)) / (m0 * Math.Pow(l21, 3)) + (km * (b4 - co * b2) * (l21 - m0) * (b3 - b1)) / (m0 * Math.Pow(l21, 3)) - (kl * (b6 - b4) * (b5 - b3)) / (l0 * Math.Pow(l22, 2)) + (kl * (b6 - b4) * (l22 - l0) * (b5 - b3)) / (l0 * Math.Pow(l22, 3)));
                Jcs.At(8 * i + 5, 8 * i + 5, -(2 * km * (l21 - m0)) / (m0 * l21) + (km * (2 * co * (co * b4 - b2) + 2 * Math.Pow(si, 2) * b4) * (co * b2 - b4)) / (2 * m0 * Math.Pow(l21, 2)) - (km * (2 * co * (co * b4 - b2) + 2 * Math.Pow(si, 2) * b4) * (b4 - co * b2)) / (2 * m0 * Math.Pow(l21, 2)) - (km * (2 * co * (co * b4 - b2) + 2 * Math.Pow(si, 2) * b4) * (co * b2 - b4) * (l21 - m0)) / (2 * m0 * Math.Pow(l21, 3)) + (km * (2 * co * (co * b4 - b2) + 2 * Math.Pow(si, 2) * b4) * (b4 - co * b2) * (l21 - m0)) / (2 * m0 * Math.Pow(l21, 3)) - (kl * (l22 - l0)) / (l0 * l22) + (kl * (b6 - b4) * (-2 * si * (si * b6 - si * b4) - 2 * co * (co * b6 - co * b4))) / (2 * l0 * Math.Pow(l22, 2)) - (kl * (b6 - b4) * (-2 * si * (si * b6 - si * b4) - 2 * co * (co * b6 - co * b4)) * (l22 - l0)) / (2 * l0 * Math.Pow(l22, 3)));
                Jcs.At(8 * i + 5, 8 * i + 6, (kl * (b6 - b4) * (b5 - b3)) / (l0 * Math.Pow(l22, 2)) - (kl * (b6 - b4) * (l22 - l0) * (b5 - b3)) / (l0 * Math.Pow(l22, 3)));
                Jcs.At(8 * i + 5, 8 * i + 7, (kl * (l22 - l0)) / (l0 * l22) + (kl * (b6 - b4) * (2 * si * (si * b6 - si * b4) + 2 * co * (co * b6 - co * b4))) / (2 * l0 * Math.Pow(l22, 2)) - (kl * (b6 - b4) * (2 * si * (si * b6 - si * b4) + 2 * co * (co * b6 - co * b4)) * (l22 - l0)) / (2 * l0 * Math.Pow(l22, 3)));

                #endregion

                /* node c */
                #region
                c1 = X[8 * i + 4]; // bx (i)
                c2 = X[8 * i + 5]; // br (i)
                c3 = X[8 * i + 6]; // cx (i)
                c4 = X[8 * i + 7]; // cr (i)
                c5 = X[8 * i + 8]; // dx (i)
                c6 = X[8 * i + 9]; // dr (i)

                l22 = Math.Sqrt(Math.Pow((c3 - c1), 2) + Math.Pow((c4 * si - c2 * si), 2) + Math.Pow((c4 * co - c2 * co), 2));  // b-c knot length
                if (l22 < l0)
                    l22 = l0;

                l33 = Math.Sqrt(Math.Pow((c5 - c3), 2) + Math.Pow((-c4 * si), 2) + Math.Pow((c6 - c4 * co), 2));    // c-d twine length
                if (l33 < m0)
                    l33 = m0;

                // Stifness matrix elements corresponding to point c - axial direction
                Jcs.At(8 * i + 6, 8 * i + 4, (kl * (l22 - l0)) / (l0 * l22) - (Math.Pow((c3 - c1), 2) * kl * (l22 - l0)) / (l0 * Math.Pow(l22, 3)) + (Math.Pow((c3 - c1), 2) * kl) / (l0 * Math.Pow(l22, 2)));
                Jcs.At(8 * i + 6, 8 * i + 5, ((c3 - c1) * kl * (-2 * si * (c4 * si - c2 * si) - 2 * co * (c4 * co - c2 * co)) * (l22 - l0)) / (2 * l0 * Math.Pow(l22, 3)) - ((c3 - c1) * kl * (-2 * si * (c4 * si - c2 * si) - 2 * co * (c4 * co - c2 * co))) / (2 * l0 * Math.Pow(l22, 2)));
                Jcs.At(8 * i + 6, 8 * i + 6, -(kl * (l22 - l0)) / (l0 * l22) + (Math.Pow((c3 - c1), 2) * kl * (l22 - l0)) / (l0 * Math.Pow(l22, 3)) - (2 * km * (l33 - m0)) / (m0 * l33) + (2 * Math.Pow((c5 - c3), 2) * km * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) - (Math.Pow((c3 - c1), 2) * kl) / (l0 * Math.Pow(l22, 2)) - (2 * Math.Pow((c5 - c3), 2) * km) / (m0 * Math.Pow(l33, 2)));
                Jcs.At(8 * i + 6, 8 * i + 7, ((c3 - c1) * kl * (2 * si * (c4 * si - c2 * si) + 2 * co * (c4 * co - c2 * co)) * (l22 - l0)) / (2 * l0 * Math.Pow(l22, 3)) - ((c5 - c3) * km * (2 * c4 * Math.Pow(si, 2) - 2 * co * (c6 - c4 * co)) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) - ((c3 - c1) * kl * (2 * si * (c4 * si - c2 * si) + 2 * co * (c4 * co - c2 * co))) / (2 * l0 * Math.Pow(l22, 2)) + ((c5 - c3) * km * (2 * c4 * Math.Pow(si, 2) - 2 * co * (c6 - c4 * co))) / (m0 * Math.Pow(l33, 2)));
                Jcs.At(8 * i + 6, 8 * i + 8, (2 * km * (l33 - m0)) / (m0 * l33) - (2 * Math.Pow((c5 - c3), 2) * km * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) + (2 * Math.Pow((c5 - c3), 2) * km) / (m0 * Math.Pow(l33, 2)));
                if (i == nx - 1 && IncludeBC)
                { }
                else
                {
                    Jcs.At(8 * i + 6, 8 * i + 9, (2 * (c5 - c3) * km * (c6 - c4 * co)) / (m0 * Math.Pow(l33, 2)) - (2 * (c5 - c3) * km * (c6 - c4 * co) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)));
                }
                // Stifness matrix elements corresponding to point c - radial direction
                Jcs.At(8 * i + 7, 8 * i + 4, ((c4 - c2) * (c3 - c1) * kl) / (l0 * Math.Pow(l22, 2)) - ((c4 - c2) * (c3 - c1) * kl * (l22 - l0)) / (l0 * Math.Pow(l22, 3)));
                Jcs.At(8 * i + 7, 8 * i + 5, (kl * (l22 - l0)) / (l0 * l22) + ((c4 - c2) * kl * (-2 * si * (c4 * si - c2 * si) - 2 * co * (c4 * co - c2 * co)) * (l22 - l0)) / (2 * l0 * Math.Pow(l22, 3)) - ((c4 - c2) * kl * (-2 * si * (c4 * si - c2 * si) - 2 * co * (c4 * co - c2 * co))) / (2 * l0 * Math.Pow(l22, 2)));
                Jcs.At(8 * i + 7, 8 * i + 6, ((c4 - c2) * (c3 - c1) * kl * (l22 - l0)) / (l0 * Math.Pow(l22, 3)) + (2 * (c5 - c3) * km * (c6 * co - c4) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) - ((c4 - c2) * (c3 - c1) * kl) / (l0 * Math.Pow(l22, 2)) - (2 * (c5 - c3) * km * (c6 * co - c4)) / (m0 * Math.Pow(l33, 2)));
                Jcs.At(8 * i + 7, 8 * i + 7, -(kl * (l22 - l0)) / (l0 * l22) + ((c4 - c2) * kl * (2 * si * (c4 * si - c2 * si) + 2 * co * (c4 * co - c2 * co)) * (l22 - l0)) / (2 * l0 * Math.Pow(l22, 3)) - (2 * km * (l33 - m0)) / (m0 * l33) - (km * (c6 * co - c4) * (2 * c4 * Math.Pow(si, 2) - 2 * co * (c6 - c4 * co)) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) - ((c4 - c2) * kl * (2 * si * (c4 * si - c2 * si) + 2 * co * (c4 * co - c2 * co))) / (2 * l0 * Math.Pow(l22, 2)) + (km * (c6 * co - c4) * (2 * c4 * Math.Pow(si, 2) - 2 * co * (c6 - c4 * co))) / (m0 * Math.Pow(l33, 2)));
                Jcs.At(8 * i + 7, 8 * i + 8, (2 * (c5 - c3) * km * (c6 * co - c4)) / (m0 * Math.Pow(l33, 2)) - (2 * (c5 - c3) * km * (c6 * co - c4) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)));
                if (i == nx - 1 && IncludeBC)
                { }
                else
                {
                    Jcs.At(8 * i + 7, 8 * i + 9, (2 * km * co * (l33 - m0)) / (m0 * l33) - (2 * km * (c6 - c4 * co) * (c6 * co - c4) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) + (2 * km * (c6 - c4 * co) * (c6 * co - c4)) / (m0 * Math.Pow(l33, 2)));
                }
                #endregion

                /*node d (first mesh)*/
                #region
                if (i == 0 && !IncludeBC)
                {
                    d3 = X[0]; // dx(0)
                    d4 = X[1]; // dr(0)
                    d5 = X[2]; // ax(1)
                    d6 = X[3]; // ar(1)

                    l44 = Math.Sqrt(Math.Pow((d5 - d3), 2) + Math.Pow((d6 - d4), 2));          // d-a knot length
                    if (l44 < l0)
                        l44 = l0;

                    // Stifness matrix elements corresponding to point d - axial direction
                    Jcs.At(0, 0, -(kl * (l44 - l0)) / (l44 * l0) + (Math.Pow((X[2] - X[0]), 2) * kl * (l44 - l0)) / (Math.Pow(l44, 3) * l0) - (Math.Pow((X[2] - X[0]), 2) * kl) / (Math.Pow(l44, 2) * l0));
                    Jcs.At(0, 1, ((X[2] - X[0]) * (X[3] - X[1]) * kl * (l44 - l0)) / (Math.Pow(l44, 3) * l0) - ((X[2] - X[0]) * (X[3] - X[1]) * kl) / (Math.Pow(l44, 2) * l0));
                    Jcs.At(0, 2, (kl * (l44 - l0)) / (l44 * l0) - (Math.Pow((X[2] - X[0]), 2) * kl * (l44 - l0)) / (Math.Pow(l44, 3) * l0) + (Math.Pow((X[2] - X[0]), 2) * kl) / (Math.Pow(l44, 2) * l0));
                    Jcs.At(0, 3, ((X[2] - X[0]) * (X[3] - X[1]) * kl) / (Math.Pow(l44, 2) * l0) - ((X[2] - X[0]) * (X[3] - X[1]) * kl * (l44 - l0)) / (Math.Pow(l44, 3) * l0));

                    // Stifness matrix elements corresponding to point d - radial direction
                    Jcs.At(1, 0, ((X[2] - X[0]) * (X[3] - X[1]) * kl * (l44 - l0)) / (Math.Pow(l44, 3) * l0) - ((X[2] - X[0]) * (X[3] - X[1]) * kl) / (Math.Pow(l44, 2) * l0));
                    Jcs.At(1, 1, -(kl * (l44 - l0)) / (l44 * l0) + (Math.Pow((X[3] - X[1]), 2) * kl * (l44 - l0)) / (Math.Pow(l44, 3) * l0) - (Math.Pow((X[3] - X[1]), 2) * kl) / (Math.Pow(l44, 2) * l0));
                    Jcs.At(1, 2, ((X[2] - X[0]) * (X[3] - X[1]) * kl) / (Math.Pow(l44, 2) * l0) - ((X[2] - X[0]) * (X[3] - X[1]) * kl * (l44 - l0)) / (Math.Pow(l44, 3) * l0));
                    Jcs.At(1, 3, (kl * (l44 - l0)) / (l44 * l0) - (Math.Pow((X[3] - X[1]), 2) * kl * (l44 - l0)) / (Math.Pow(l44, 3) * l0) + (Math.Pow((X[3] - X[1]), 2) * kl) / (Math.Pow(l44, 2) * l0));

                }
                #endregion

                /*node d (last mesh)*/
                #region
                if (i == nx - 1)
                {
                    d1 = X[8 * i + 6]; // cx(last)
                    d2 = X[8 * i + 7]; // cr(last)
                    d3 = X[8 * i + 8]; // dx(last)
                    d4 = X[8 * i + 9]; // dr(last)

                    l33 = Math.Sqrt(Math.Pow((d3 - d1), 2) + Math.Pow((-d2 * si), 2) + Math.Pow((d4 - d2 * co), 2));                // c-d knot length
                    if (l33 < m0)
                        l33 = m0;

                    // Stifness matrix elements corresponding to point d - axial direction
                    Jcs.At(8 * i + 8, 8 * i + 6, (2 * km * (l33 - m0)) / (m0 * l33) - (Math.Pow((d3 - d1), 2) * km * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) + ((d1 - d3) * (d3 - d1) * km * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) + (Math.Pow((d3 - d1), 2) * km) / (m0 * Math.Pow(l33, 2)) - ((d1 - d3) * (d3 - d1) * km) / (m0 * Math.Pow(l33, 2)));
                    Jcs.At(8 * i + 8, 8 * i + 7, ((d3 - d1) * km * (2 * d2 * Math.Pow(si, 2) - 2 * co * (d4 - d2 * co)) * (l33 - m0)) / (2 * m0 * Math.Pow(l33, 3)) - ((d1 - d3) * km * (2 * d2 * Math.Pow(si, 2) - 2 * co * (d4 - d2 * co)) * (l33 - m0)) / (2 * m0 * Math.Pow(l33, 3)) - ((d3 - d1) * km * (2 * d2 * Math.Pow(si, 2) - 2 * co * (d4 - d2 * co))) / (2 * m0 * Math.Pow(l33, 2)) + ((d1 - d3) * km * (2 * d2 * Math.Pow(si, 2) - 2 * co * (d4 - d2 * co))) / (2 * m0 * Math.Pow(l33, 2)));
                    Jcs.At(8 * i + 8, 8 * i + 8, -(2 * km * (l33 - m0)) / (m0 * l33) + (Math.Pow((d3 - d1), 2) * km * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) - ((d1 - d3) * (d3 - d1) * km * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) - (Math.Pow((d3 - d1), 2) * km) / (m0 * Math.Pow(l33, 2)) + ((d1 - d3) * (d3 - d1) * km) / (m0 * Math.Pow(l33, 2)));
                    if (!IncludeBC)
                    {
                        Jcs.At(8 * i + 8, 8 * i + 9, ((d3 - d1) * km * (d4 - d2 * co) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) - ((d1 - d3) * km * (d4 - d2 * co) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) - ((d3 - d1) * km * (d4 - d2 * co)) / (m0 * Math.Pow(l33, 2)) + ((d1 - d3) * km * (d4 - d2 * co)) / (m0 * Math.Pow(l33, 2)));
                    }
                    // Stifness matrix elements corresponding to point d - radial direction
                    if (!IncludeBC)
                    {
                        Jcs.At(8 * i + 9, 8 * i + 6, ((d3 - d1) * km * (d2 * co - d4) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) - ((d3 - d1) * km * (d4 - d2 * co) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) - ((d3 - d1) * km * (d2 * co - d4)) / (m0 * Math.Pow(l33, 2)) + ((d3 - d1) * km * (d4 - d2 * co)) / (m0 * Math.Pow(l33, 2)));
                        Jcs.At(8 * i + 9, 8 * i + 7, (2 * km * co * (l33 - m0)) / (m0 * l33) - (km * (d2 * co - d4) * (2 * d2 * Math.Pow(si, 2) - 2 * co * (d4 - d2 * co)) * (l33 - m0)) / (2 * m0 * Math.Pow(l33, 3)) + (km * (d4 - d2 * co) * (2 * d2 * Math.Pow(si, 2) - 2 * co * (d4 - d2 * co)) * (l33 - m0)) / (2 * m0 * Math.Pow(l33, 3)) + (km * (d2 * co - d4) * (2 * d2 * Math.Pow(si, 2) - 2 * co * (d4 - d2 * co))) / (2 * m0 * Math.Pow(l33, 2)) - (km * (d4 - d2 * co) * (2 * d2 * Math.Pow(si, 2) - 2 * co * (d4 - d2 * co))) / (2 * m0 * Math.Pow(l33, 2)));
                        Jcs.At(8 * i + 9, 8 * i + 8, -((d3 - d1) * km * (d2 * co - d4) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) + ((d3 - d1) * km * (d4 - d2 * co) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) + ((d3 - d1) * km * (d2 * co - d4)) / (m0 * Math.Pow(l33, 2)) - ((d3 - d1) * km * (d4 - d2 * co)) / (m0 * Math.Pow(l33, 2)));
                        Jcs.At(8 * i + 9, 8 * i + 9, -(2 * km * (l33 - m0)) / (m0 * l33) - (km * (d4 - d2 * co) * (d2 * co - d4) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) + (km * Math.Pow((d4 - d2 * co), 2) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) + (km * (d4 - d2 * co) * (d2 * co - d4)) / (m0 * Math.Pow(l33, 2)) - (km * Math.Pow((d4 - d2 * co), 2)) / (m0 * Math.Pow(l33, 2)));
                    }

                }
                #endregion

                /* node d (intermediate meshes)*/
                #region
                else
                {
                    d1 = X[8 * i + 6]; // cx (i)
                    d2 = X[8 * i + 7]; // cr (i)
                    d3 = X[8 * i + 8]; // dx (i)
                    d4 = X[8 * i + 9]; // dr (i)
                    d5 = X[8 * i + 10]; // ax (i+1)
                    d6 = X[8 * i + 11]; // ar (i+1)

                    l33 = Math.Sqrt(Math.Pow((d3 - d1), 2) + Math.Pow((-d2 * si), 2) + Math.Pow((d4 - d2 * co), 2));        // c-d twine length
                    if (l33 < m0)
                        l33 = m0;

                    l44 = Math.Sqrt(Math.Pow((d5 - d3), 2) + Math.Pow((d6 - d4), 2));                                   // d-a knot length
                    if (l44 < l0)
                        l44 = l0;

                    // Stifness matrix elements corresponding to point d - axial direction
                    Jcs.At(8 * i + 8, 8 * i + 6, (2 * km * (l33 - m0)) / (m0 * l33) - (Math.Pow((d3 - d1), 2) * km * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) + ((d1 - d3) * (d3 - d1) * km * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) + (Math.Pow((d3 - d1), 2) * km) / (m0 * Math.Pow(l33, 2)) - ((d1 - d3) * (d3 - d1) * km) / (m0 * Math.Pow(l33, 2)));
                    Jcs.At(8 * i + 8, 8 * i + 7, ((d3 - d1) * km * (2 * d2 * Math.Pow(si, 2) - 2 * co * (d4 - d2 * co)) * (l33 - m0)) / (2 * m0 * Math.Pow(l33, 3)) - ((d1 - d3) * km * (2 * d2 * Math.Pow(si, 2) - 2 * co * (d4 - d2 * co)) * (l33 - m0)) / (2 * m0 * Math.Pow(l33, 3)) - ((d3 - d1) * km * (2 * d2 * Math.Pow(si, 2) - 2 * co * (d4 - d2 * co))) / (2 * m0 * Math.Pow(l33, 2)) + ((d1 - d3) * km * (2 * d2 * Math.Pow(si, 2) - 2 * co * (d4 - d2 * co))) / (2 * m0 * Math.Pow(l33, 2)));
                    Jcs.At(8 * i + 8, 8 * i + 8, -(2 * km * (l33 - m0)) / (m0 * l33) + (Math.Pow((d3 - d1), 2) * km * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) - ((d1 - d3) * (d3 - d1) * km * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) - (Math.Pow((d3 - d1), 2) * km) / (m0 * Math.Pow(l33, 2)) + ((d1 - d3) * (d3 - d1) * km) / (m0 * Math.Pow(l33, 2)) - (kl * (l44 - l0)) / (l44 * l0) + (Math.Pow((d5 - d3), 2) * kl * (l44 - l0)) / (Math.Pow(l44, 3) * l0) - (Math.Pow((d5 - d3), 2) * kl) / (Math.Pow(l44, 2) * l0));
                    Jcs.At(8 * i + 8, 8 * i + 9, ((d3 - d1) * km * (d4 - d2 * co) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) - ((d1 - d3) * km * (d4 - d2 * co) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) - ((d3 - d1) * km * (d4 - d2 * co)) / (m0 * Math.Pow(l33, 2)) + ((d1 - d3) * km * (d4 - d2 * co)) / (m0 * Math.Pow(l33, 2)) + ((d5 - d3) * (d6 - d4) * kl * (l44 - l0)) / (Math.Pow(l44, 3) * l0) - ((d5 - d3) * (d6 - d4) * kl) / (Math.Pow(l44, 2) * l0));
                    Jcs.At(8 * i + 8, 8 * i + 10, (kl * (l44 - l0)) / (l44 * l0) - (Math.Pow((d5 - d3), 2) * kl * (l44 - l0)) / (Math.Pow(l44, 3) * l0) + (Math.Pow((d5 - d3), 2) * kl) / (Math.Pow(l44, 2) * l0));
                    Jcs.At(8 * i + 8, 8 * i + 11, ((d5 - d3) * (d6 - d4) * kl) / (Math.Pow(l44, 2) * l0) - ((d5 - d3) * (d6 - d4) * kl * (l44 - l0)) / (Math.Pow(l44, 3) * l0));

                    // Stifness matrix elements corresponding to point d - axial direction
                    Jcs.At(8 * i + 9, 8 * i + 6, -((d3 - d1) * (d4 - co * d2) * km * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) + ((d3 - d1) * (co * d2 - d4) * km * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) + ((d3 - d1) * (d4 - co * d2) * km) / (m0 * Math.Pow(l33, 2)) - ((d3 - d1) * (co * d2 - d4) * km) / (m0 * Math.Pow(l33, 2)));
                    Jcs.At(8 * i + 9, 8 * i + 7, (2 * co * km * (l33 - m0)) / (m0 * l33) + ((d4 - co * d2) * km * (2 * d2 * Math.Pow(si, 2) - 2 * co * (d4 - co * d2)) * (l33 - m0)) / (2 * m0 * Math.Pow(l33, 3)) - ((co * d2 - d4) * km * (2 * d2 * Math.Pow(si, 2) - 2 * co * (d4 - co * d2)) * (l33 - m0)) / (2 * m0 * Math.Pow(l33, 3)) - ((d4 - co * d2) * km * (2 * d2 * Math.Pow(si, 2) - 2 * co * (d4 - co * d2))) / (2 * m0 * Math.Pow(l33, 2)) + ((co * d2 - d4) * km * (2 * d2 * Math.Pow(si, 2) - 2 * co * (d4 - co * d2))) / (2 * m0 * Math.Pow(l33, 2)));
                    Jcs.At(8 * i + 9, 8 * i + 8, -((d3 - d1) * km * (d2 * co - d4) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) + ((d3 - d1) * km * (d4 - d2 * co) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) + ((d3 - d1) * km * (d2 * co - d4)) / (m0 * Math.Pow(l33, 2)) - ((d3 - d1) * km * (d4 - d2 * co)) / (m0 * Math.Pow(l33, 2)) + ((d5 - d3) * (d6 - d4) * kl * (l44 - l0)) / (Math.Pow(l44, 3) * l0) - ((d5 - d3) * (d6 - d4) * kl) / (Math.Pow(l44, 2) * l0));
                    Jcs.At(8 * i + 9, 8 * i + 9, -(2 * km * (l33 - m0)) / (m0 * l33) - (km * (d4 - d2 * co) * (d2 * co - d4) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) + (km * Math.Pow((d4 - d2 * co), 2) * (l33 - m0)) / (m0 * Math.Pow(l33, 3)) + (km * (d4 - d2 * co) * (d2 * co - d4)) / (m0 * Math.Pow(l33, 2)) - (km * Math.Pow((d4 - d2 * co), 2)) / (m0 * Math.Pow(l33, 2)) - (kl * (l44 - l0)) / (l44 * l0) + (Math.Pow((d6 - d4), 2) * kl * (l44 - l0)) / (Math.Pow(l44, 3) * l0) - (Math.Pow((d6 - d4), 2) * kl) / (Math.Pow(l44, 2) * l0));
                    Jcs.At(8 * i + 9, 8 * i + 10, ((d5 - d3) * (d6 - d4) * kl) / (Math.Pow(l44, 2) * l0) - ((d5 - d3) * (d6 - d4) * kl * (l44 - l0)) / (Math.Pow(l44, 3) * l0));
                    Jcs.At(8 * i + 9, 8 * i + 11, (kl * (l44 - l0)) / (l44 * l0) - (Math.Pow((d6 - d4), 2) * kl * (l44 - l0)) / (Math.Pow(l44, 3) * l0) + (Math.Pow((d6 - d4), 2) * kl) / (Math.Pow(l44, 2) * l0));

                    #endregion
                }
            }
        }

        private void SetCatchJacobian(bool IncludeBC)
        {
            /* intermediate meshes */
            #region
            for (int i = ncp + 1; i < dof / 2; i++)
            {
                // Stffness correpsonding to axial force component 
                Jcs.At(2 * i - 2, 2 * i - 3, P * theta * X[2 * i - 3]);
                if (i == dof / 2 - 1 && IncludeBC)
                { }
                else
                {
                    Jcs.At(2 * i - 2, 2 * i + 1, -P * theta * X[2 * i + 1]);
                }
                // Stffness correpsonding to radial force component 
                Jcs.At(2 * i - 1, 2 * i - 4, -0.5 * P * theta * (X[2 * i - 3] + X[2 * i - 1]) * (X[2 * i - 2] - X[2 * i - 4]) / Math.Abs(X[2 * i - 2] - X[2 * i - 4]));
                Jcs.At(2 * i - 1, 2 * i - 3, 0.5 * P * theta * Math.Abs(X[2 * i - 2] - X[2 * i - 4]));
                Jcs.At(2 * i - 1, 2 * i - 2, 0.5 * P * theta * (X[2 * i - 3] + X[2 * i - 1]) * (X[2 * i - 2] - X[2 * i - 4]) / Math.Abs(X[2 * i - 2] - X[2 * i - 4]) -
                                               0.5 * P * theta * (X[2 * i] - X[2 * i - 2]) * (X[2 * i + 1] + X[2 * i - 1]) / Math.Abs(X[2 * i] - X[2 * i - 2]));
                Jcs.At(2 * i - 1, 2 * i - 1, 0.5 * P * theta * Math.Abs(X[2 * i] - X[2 * i - 2]) + 0.5 * P * theta * Math.Abs(X[2 * i - 2] - X[2 * i - 4]));
                Jcs.At(2 * i - 1, 2 * i - 0, 0.5 * P * theta * (X[2 * i] - X[2 * i - 2]) * (X[2 * i + 1] + X[2 * i - 1]) / Math.Abs(X[2 * i] - X[2 * i - 2]));
                if (i == dof / 2 - 1 && IncludeBC)
                {
                }
                else
                {
                    Jcs.At(2 * i - 1, 2 * i + 1, 0.5 * P * theta * Math.Abs(X[2 * i] - X[2 * i - 2]));
                }
            }
            #endregion

            /*first mesh*/
            #region
            if (ncp == 1 && IncludeBC)
            { }
            else
            {
                // Stffness correpsonding to axial force component
                Jcs.At(2 * ncp - 2, 2 * ncp - 1, P * theta * X[2 * ncp - 1]);
                Jcs.At(2 * ncp - 2, 2 * ncp + 1, -P * theta * X[2 * ncp + 1]);
                // Stffness correpsonding to radial force component 
                Jcs.At(2 * ncp - 1, 2 * ncp - 2, -0.5 * P * theta * (X[2 * ncp] - X[2 * ncp - 2]) * (X[2 * ncp + 1] + X[2 * ncp - 1]) / Math.Abs(X[2 * ncp] - X[2 * ncp - 2]));
                Jcs.At(2 * ncp - 1, 2 * ncp - 1, 0.5 * P * theta * Math.Abs(X[2 * ncp] - X[2 * ncp - 2]));
                Jcs.At(2 * ncp - 1, 2 * ncp - 0, 0.5 * P * theta * (X[2 * ncp] - X[2 * ncp - 2]) * (X[2 * ncp + 1] + X[2 * ncp - 1]) / Math.Abs(X[2 * ncp] - X[2 * ncp - 2]));
                Jcs.At(2 * ncp - 1, 2 * ncp + 1, 0.5 * P * theta * Math.Abs(X[2 * ncp] - X[2 * ncp - 2]));
            }
            #endregion

            /*last mesh*/
            #region
            // Stffness correpsonding to axial force component 
            Jcs.At(dof - 2, dof - 3, P * theta * X[dof - 3]);
            if (!IncludeBC)
            {
                Jcs.At(dof - 2, dof - 1, -P * theta * X[dof - 1]);
            }
            // Stffness correpsonding to radial force component 
            if (!IncludeBC)
            {
                Jcs.At(dof - 1, dof - 4, -0.5 * P * theta * (X[dof - 2] - X[dof - 4]) * (X[dof - 1] + X[dof - 3]) / Math.Abs(X[dof - 2] - X[dof - 4]));
                Jcs.At(dof - 1, dof - 3, 0.5 * P * theta * Math.Abs(X[dof - 2] - X[dof - 4]));
                Jcs.At(dof - 1, dof - 2, 0.5 * P * theta * (X[dof - 2] - X[dof - 4]) * (X[dof - 1] + X[dof - 3]) / Math.Abs(X[dof - 2] - X[dof - 4]));
                Jcs.At(dof - 1, dof - 1, 0.5 * P * theta * Math.Abs(X[dof - 2] - X[dof - 4]));
            }
            #endregion
        }

        public void UpdateTotalJacobian(double kDiag, bool IncludeBC)
        {
            Jcs = new CoordinateStorage<double>(dof, dof, 4 * 4 + (dof - 4) * 6);
            SetTwineJacobian(IncludeBC);
            SetCatchJacobian(IncludeBC);
            Jcs.At(0, 0, -1);
            Jcs.At(1, 1, -1);
            Jcs.At(dof - 1, dof - 1, -1);

            if (kDiag != 0)
            {
                for (int i = 0; i < dof; i++)
                {
                    Jcs.At(i, i, kDiag);
                }
            }

            J = Converter.ToCompressedColumnStorage(Jcs);
        }

        /* results*/

        public double Length()
        {
            return X[dof - 2];
        }

        public double MaxRadius()
        {
            double rMax = r0;                           // maximum radius
            for (int i = 1; i < dof / 2; i++)
            {
                if (rMax < X[2 * i - 1])
                {
                    rMax = X[2 * i - 1];
                }
            }
            return rMax;
        }

        public double EntranceDrag()
        {
            UpdateTotalForces(IncludeBC: false);
            return nr * F[0];                                   // resultant entrance reaction
        }

        public double CatchDrag()
        {
            double rCatch = X[2 * ncp - 1];      // radius at the beginning of the catch
            return pi * Math.Pow(rCatch, 2) * P; // resultant catch drag
        }

        public double CatchThickness()
        {
            return X[dof - 2] - X[2 * ncp - 2];
        }

        public double CatchVolume()
        {
            double V = 0;
            for (int i = ncp; i < dof / 2; i++)
            {                                           
                V = V + pi * (
                              X[2 * i] - X[2 * i - 2]) * (Math.Pow(X[2 * i + 1], 2) + 
                              X[2 * i + 1] * X[2 * i - 1] + Math.Pow(X[2 * i - 1], 2)
                              ) / 3; // volume of a cut cone between 2 nodes
            }
            return V;
        }

        public double CatchSurface()
        {
            double S = 0;
            for (int i = ncp; i < dof / 2; i++)         
            {
                S = S + pi * (X[2 * i - 1] + X[2 * i + 1]) * Math.Pow(Math.Pow(X[2 * i] - X[2 * i - 2], 2) + 
                              Math.Pow(X[2 * i + 1] - X[2 * i - 1], 2), 0.5); // surface area of a cut cone between 2 nodes
            }
            return S;
        }

        /* printing*/

        public void PrintResults()
        {
            Console.WriteLine("\nResults:");
            Console.WriteLine("Total codend length:             {0,10:N3} [m]", X[dof - 2]);
            Console.WriteLine("Maximum codend radius:           {0,10:N3} [m]", MaxRadius());
            Console.WriteLine("Catch extent:                    {0,10:N3} [m]", CatchThickness());
            Console.WriteLine("Surface in contact with catch:   {0,10:N3} [m^2]", CatchSurface());
            Console.WriteLine("Catch volume:                    {0,10:N3} [m^3]", CatchVolume());
            Console.WriteLine("Resultant entrance reaction:     {0,10:N0} [N]", EntranceDrag());
            Console.WriteLine("Total catch drag force:          {0,10:N0} [N]\n", CatchDrag());
        }

        private void PrintVector(double[] Array)
        {
            Console.WriteLine();

            for (int i = 0; i < Array.Length; i++)
            {
                Console.WriteLine(String.Format("{0,-15:E4}", Array[i]));
            }
        }

        private void PrintMatrix(double[,] Array, int[] Rows, int[] Cols)
        {
            Console.WriteLine();
            for (int i = 0; i < Rows.Length; i++)
            {
                for (int j = 0; j < Cols.Length; j++)
                {
                    Console.Write(String.Format("{0,-15:E4}", Array[Rows[i], Cols[j]]));
                }
                Console.WriteLine();
            }
        }
      
    } 
}
