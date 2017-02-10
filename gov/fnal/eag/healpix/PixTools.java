package gov.fnal.eag.healpix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Vector;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.lang.Number;

import javax.vecmath.Vector3d;

/**
 * 
 *  contains methods translated from HEALPix Fortran90
 *  with increased map resolution in comparison to original Fortran code.
 * 
 * @author N Kuropatkin
 * 
 * Created on Mar 10, 2005
 * Modified on December 18 2007
 * Corrected arithmetic and mistyping April 20 2008
 * 
 * @author Mark Taylor made modifications to make the class thread safe 11-Jan-2012
 * 
 * <p>
 *  All methods are thread safe.  This class can be used as a singleton,
 *  the singleton instance being available from the {@link #getInstance} method.
 *  For compatibility with previous versions however it is possible to
 *  construct new instances using the default constructor.
 *</p>
 *  
 */
public class PixTools {

	private static final double twothird = 2. / 3.;

	private static final double PI = Math.PI;

	private static final double TWOPI = 2. * PI;



	private static final double HALFPI = PI / 2.0;

    private static final int ns_max = 1048576; // 2^20

//
    private static final int xmax = 4096;
//
    private static final int pixmax = 262144;

//
    private static final int xmid = 512;   
    
	private static final long[] x2pix = new long[xmax+1];

	private static final long[] y2pix = new long[xmax+1];

	private static final long[] pix2x = new long[pixmax+1];

	private static final long[] pix2y = new long[pixmax+1];



	static {
		mk_xy2pix();
	}

        /** Singleton instance. */
	private static final PixTools pixTools = new PixTools();

	/**
	 * default constructor
	 * 
	 *  
	 */
	public PixTools() {
	}

	/**
	 * finds pixels having a colatitude (measured from North pole) : 
	 * theta1 < colatitude < theta2 with 0 <= theta1 < theta2 <= Pi 
	 * if theta2 < theta1
	 * then pixels with 0 <= colatitude < theta2 or theta1 < colatitude < Pi are
	 * returned
	 * 
	 * @param nside 
	 *            long the map resolution parameter
	 * @param theta1 
	 *            lower edge of the colatitude
	 * @param theta2 
	 *            upper edge of the colatitude
	 * @param nest 
	 *            long if = 1 result is in NESTED scheme
	 * @return  ArrayList of  pixel numbers (long)
	 * @throws Exception 
	 * @throws IllegalArgumentException
	 */
	public ArrayList query_strip(long nside, double theta1, double theta2,
			long nest) throws Exception {
		ArrayList res = new ArrayList();
		ArrayList listir = new ArrayList();
		long npix, nstrip;
		long iz,  irmin, irmax;
        int is;
		double phi0, dphi;
		double[] colrange = new double[4];
		boolean nest_flag = false;
		String SID = " QUERY_STRIP";
		/* ---------------------------------------- */
		npix = Nside2Npix(nside);
		if (nest == 1)
			nest_flag = true;
		if (npix < 0) {
			throw new IllegalArgumentException(SID + " Nside should be power of 2");
		}
		if ((theta1 < 0.0 || theta1 > PI) || (theta2 < 0.0 || theta2 > PI)) {
			throw new IllegalArgumentException(SID + " Illegal value of theta1, theta2");
		}
		if (theta1 <= theta2) {
			nstrip = 1;
			colrange[0] = theta1;
			colrange[1] = theta2;
		} else {
			nstrip = 2;
			colrange[0] = 0.0;
			colrange[1] = theta2;
			colrange[2] = theta1;
			colrange[3] = PI;
		}
		/* loops on strips */
		for (is = 0; is < nstrip; is++) {
			irmin = RingNum(nside, Math.cos(colrange[2 * is]));
			irmax = RingNum(nside, Math.cos(colrange[2 * is + 1]));
			/* loop on ring number */
			for (iz = irmin; iz <= irmax; iz++) {
				phi0 = 0.;
				dphi = PI;
				listir = InRing(nside, iz, phi0, dphi, nest_flag);
				res.addAll(listir);
			}
		}
		return res;
	}

	/**
	 * finds pixels that lay within a CONVEX polygon defined by its vertex on
	 * sphere
	 * 
	 * @param nside 
	 *            the map resolution
	 * @param vlist 
	 *            ArrayList of vectors defining the polygon vertices
	 * @param nest 
	 *            if set to 1 use NESTED scheme
	 * @param inclusive 
	 *            if set 1 returns all pixels crossed by polygon boundaries
	 * @return  ArrayList of pixels
	 * 
	 * algorithm: the polygon is divided into triangles vertex 0 belongs to all
	 * triangles
	 * @throws Exception 
	 * @throws IllegalArgumentException
	 */
	public ArrayList query_polygon(long nside, ArrayList vlist, long nest,
			long inclusive) throws Exception {
		ArrayList res = new ArrayList();
		int nv = vlist.size();
		Vector3d vp0, vp1, vp2;
		Vector3d vo;
		ArrayList vvlist = new ArrayList();
//		double surface, fsky;
		double hand;
		double[] ss = new double[nv];
//		int n_in_trg, ilist, ntl;
        long npix;
		int ix = 0;
		
		int n_remain, np, nm, nlow;
		String SID = "QUERY_POLYGON";

		//		System.out.println("Start polygon");
		for (int k = 0; k < nv; k++)
			ss[k] = 0.;
		/* -------------------------------------- */
		n_remain = nv;
		if (n_remain < 3) {
			throw new IllegalArgumentException(SID + " Number of vertices should be >= 3");
		}
		/*---------------------------------------------------------------- */
		/* Check that the poligon is convex or has only one concave vertex */
		/*---------------------------------------------------------------- */
		int i0 = 0;
		int i2 = 0;
		if (n_remain > 3) { // a triangle is always convex
			for (int i1 = 1; i1 <= n_remain - 1; i1++) { // in [0,n_remain-1]
				i0 = (int) BitManipulation.MODULO(i1 - 1, n_remain);
				i2 = (int) BitManipulation.MODULO(i1 + 1, n_remain);
				vp0 = (Vector3d) vlist.get(i0); // select vertices by 3
												// neighbour
				vp1 = (Vector3d) vlist.get(i1);
				vp2 = (Vector3d) vlist.get(i2);
				// computes handedness (v0 x v2) . v1 for each vertex v1
				vo = new Vector3d(crossProduct(vp0, vp2));
				hand = dotProduct(vo, vp1);
				if (hand >= 0.) {
					ss[i1] = 1.0;
				} else {
					ss[i1] = -1.0;
				}

			}
			np = 0; // number of vert. with positive handedness
			for (int i = 0; i < nv; i++) {
				if (ss[i] > 0.)
					np++;
			}
			nm = n_remain - np;

			nlow = Math.min(np, nm);

			if (nlow != 0) {
				if (nlow == 1) { // only one concave vertex
					if (np == 1) { // ix index of the vertex in the list
						for (int k = 0; k < nv - 1; k++) {
							if (Math.abs(ss[k] - 1.0) <= 1.e-12) {
								ix = k;
								break;
							}
						}
					} else {
						for (int k = 0; k < nv - 1; k++) {
							if (Math.abs(ss[k] + 1.0) <= 1.e-12) {
								ix = k;
								break;
							}
						}
					}

					// rotate pixel list to put that vertex in #0
					int n_rot = vlist.size() - ix;
					int ilast = vlist.size() - 1;
					for (int k = 0; k < n_rot; k++) {
						Vector3d temp = new Vector3d((Vector3d) vlist
								.get(ilast));
						vlist.remove(ilast);
						vlist.add(0,  temp);
					}
				}
				if (nlow > 1) { // more than 1concave vertex
					System.out
							.println(" The polygon has more than one concave vertex");
					System.out.println(" The result is unpredictable");
				}
			}
		}
		/* fill the polygon, one triangle at a time */
		npix = (long) Nside2Npix(nside);
		while (n_remain >= 3) {
			vp0 = (Vector3d) vlist.get(0);
			vp1 = (Vector3d) vlist.get(n_remain - 2);
			vp2 = (Vector3d) vlist.get(n_remain - 1);

			/* find pixels within the triangle */
			ArrayList templist = new ArrayList();
			templist = query_triangle(nside, vp0, vp1, vp2, nest, inclusive);

			vvlist.addAll(templist);
			n_remain--;
		}
		/* make final pixel list */
		npix = vvlist.size();
		long[] pixels = new long[(int)npix];
		for (int i = 0; i < npix; i++) {
			pixels[i] = ((Long) vvlist.get(i)).longValue();
		}
		Arrays.sort(pixels);
		int k = 0;
		res.add(k, new Long(pixels[0]));
		for (int i = 1; i < pixels.length; i++) {
			if (pixels[i] > pixels[i - 1]) {
				k++;
				res.add(k, new Long(pixels[i]));
			}
		}

		return res;
	}

	/**
	 * generates a list of pixels that lay inside a triangle defined by
	 * the three vertex vectors
	 * 
	 * @param nside 
	 *            long map resolution parameter
	 * @param v1 
	 *            Vector3d defines one vertex of the triangle
	 * @param v2 
	 *            Vector3d another vertex
	 * @param v3 
	 *            Vector3d yet another one
	 * @param nest 
	 *            long 0 (default) RING numbering scheme, if set to 1 the NESTED
	 *            scheme will be used.
	 * @param inclusive
	 *            long 0 (default) only pixels whose centers are inside the
	 *            triangle will be listed, if set to 1 all pixels overlaping the
	 *            triangle will be listed
	 * @return ArrayList with pixel numbers
	 * @throws Exception 
	 * @throws IllegalArgumentException
	 */
	public ArrayList query_triangle(long nside, Vector3d v1, Vector3d v2,
			Vector3d v3, long nest, long inclusive) throws Exception {
		ArrayList res;
		res = new ArrayList();
		ArrayList listir;
		long npix, iz, irmin, irmax, n12, n123a, n123b, ndom = 0;
		boolean test1, test2, test3;
		boolean test1a, test1b, test2a, test2b, test3a, test3b;
		double dth1, dth2, determ, sdet;
		double zmax, zmin, z1max, z1min, z2max, z2min, z3max, z3min;
		double z, tgth, st, offset, sin_off;
		double phi_pos, phi_neg;
		Vector3d[] vv = new Vector3d[3];
		Vector3d[] vo = new Vector3d[3];
		double[] sprod = new double[3];
		double[] sto = new double[3];
		double[] phi0i = new double[3];
		double[] tgthi = new double[3];
		double[] dc = new double[3];
		double[][] dom = new double[3][2];
		double[] dom12 = new double[4];
		double[] dom123a = new double[4];
		double[] dom123b = new double[4];
		double[] alldom = new double[6];
		double a_i, b_i, phi0, dphiring;
		long idom;
//		long nir, ip, status;
		boolean do_inclusive = false;
		boolean do_nest = false;
		String SID = "QUERY_TRIANGLE";
		long nsidesq = nside * nside;
		/*                                       */

		//		System.out.println("in query_triangle");
		npix = Nside2Npix(nside);
		if (npix < 0) {
			throw new IllegalArgumentException(SID + " Nside should be power of 2 >0 and < "+ns_max);
		}
		if (inclusive == 1)
			do_inclusive = true;
		if (nest == 1)
			do_nest = true;
		vv[0] = new Vector3d(v1);
		vv[0].normalize();
		vv[1] = new Vector3d(v2);
		vv[1].normalize();
		vv[2] = new Vector3d(v3);
		vv[2].normalize();

		/*                                  */
		dth1 = 1.0 / (3.0 * nsidesq);
		dth2 = 2.0 / (3.0 * nside);
		/*
		 * determ = (v1 X v2) . v3 determines the left ( <0) or right (>0)
		 * handedness of the triangle
		 */
		Vector3d vt = new Vector3d(0., 0., 0.);
		vt = crossProduct(vv[0], vv[1]);
		determ = dotProduct(vt, vv[2]);

		if (Math.abs(determ) < 1.0e-20) {
			throw new HealpixException(
					SID
							+ ": the triangle is degenerated - query cannot be performed");
		}
		if (determ >= 0.) { // The sign of determinant
			sdet = 1.0;
		} else {
			sdet = -1.0;
		}

		sprod[0] = dotProduct(vv[1], vv[2]);
		sprod[1] = dotProduct(vv[2], vv[0]);
		sprod[2] = dotProduct(vv[0], vv[1]);
		/* vector orthogonal to the great circle containing the vertex doublet */

		vo[0] = crossProduct(vv[1], vv[2]);
		vo[1] = crossProduct(vv[2], vv[0]);
		vo[2] = crossProduct(vv[0], vv[1]);
		vo[0].normalize();
		vo[1].normalize();
		vo[2].normalize();

		/* test presence of poles in the triangle */
		zmax = -1.0;
		zmin = 1.0;
		test1 = (vo[0].z * sdet >= 0.0); // north pole in hemisphere defined by
		// 2-3
		test2 = (vo[1].z * sdet >= 0.0); // north pole in the hemisphere defined
		// by 1-2
		test3 = (vo[2].z * sdet >= 0.0); // north pole in hemisphere defined by
		// 1-3
		if (test1 && test2 && test3)
			zmax = 1.0; // north pole in the triangle
		if ((!test1) && (!test2) && (!test3))
			zmin = -1.0; // south pole in the triangle
		/* look for northenest and southernest points in the triangle */
		test1a = ((vv[2].z - sprod[0] * vv[1].z) >= 0.0); // segment 2-3
		test1b = ((vv[1].z - sprod[0] * vv[2].z) >= 0.0);
		test2a = ((vv[2].z - sprod[1] * vv[0].z) >= 0.0); // segment 1-3
		test2b = ((vv[0].z - sprod[1] * vv[2].z) >= 0.0);
		test3a = ((vv[1].z - sprod[2] * vv[0].z) >= 0.0); // segment 1-2
		test3b = ((vv[0].z - sprod[2] * vv[1].z) >= 0.0);

		/* sin of theta for orthogonal vector */
		for (int i = 0; i < 3; i++) {
			sto[i] = Math.sqrt((1.0 - vo[i].z) * (1.0 + vo[i].z));
		}
		/*
		 * for each segment ( side of the triangle ) the extrema are either -
		 * -the 2 vertices - one of the vertices and a point within the segment
		 */
		// segment 2-3
		z1max = vv[1].z;
		z1min = vv[2].z;
		//		if (test1a == test1b) {
		//			zz = sto[0];
		//			if (vv[1].z + vv[2].z >= 0.0) {
		//				z1max = zz;
		//			} else {
		//				z1min = -zz;
		//			}
		//		}
		// segment 1-3
		z2max = vv[2].z;
		z2min = vv[0].z;
		//		if (test2a == test2b) {
		//			zz = sto[1];
		//			if (vv[0].z + vv[2].z >= 0.0) {
		//				z2max = zz;
		//			} else {
		//				z2min = -zz;
		//			}
		//		}
		// segment 1-2
		z3max = vv[0].z;
		z3min = vv[1].z;
		//		if (test3a == test3b) {
		//			zz = sto[2];
		//			if (vv[0].z + vv[1].z >= 0.0) {
		//				z3max = zz;
		//			} else {
		//				z3min = -zz;
		//			}
		//		}

		zmax = Math.max(Math.max(z1max, z2max), Math.max(z3max, zmax));
		zmin = Math.min(Math.min(z1min, z2min), Math.min(z3min, zmin));
		/*
		 * if we are inclusive, move upper point up, and lower point down, by a
		 * half pixel size
		 */
		offset = 0.0;
		sin_off = 0.0;
		if (do_inclusive) {
			offset = PI / (4.0 * nside); // half pixel size
			sin_off = Math.sin(offset);
			zmax = Math.min(1.0, Math.cos(Math.acos(zmax) - offset));
			zmin = Math.max(-1.0, Math.cos(Math.acos(zmin) + offset));
		}

		irmin = RingNum(nside, zmax);
		irmax = RingNum(nside, zmin);

		//		System.out.println("irmin = " + irmin + " irmax =" + irmax);

		/* loop on the rings */
		for (int i = 0; i < 3; i++) {
			tgthi[i] = -1.0e30 * vo[i].z;
			phi0i[i] = 0.0;
		}
		for (int j = 0; j < 3; j++) {
			if (sto[j] > 1.0e-10) {
				tgthi[j] = -vo[j].z / sto[j]; // - cotan(theta_orth)

				phi0i[j] = Math.atan2(vo[j].y, vo[j].x); // Should make it 0-2pi
														 // ?
				/* Bring the phi0i to the [0,2pi] domain if need */

				if (phi0i[j] < 0.) {
					phi0i[j] = BitManipulation.MODULO(
							(Math.atan2(vo[j].y, vo[j].x) + TWOPI), TWOPI); //  [0-2pi]
				}

			}
		}
		/*
		 * the triangle boundaries are geodesics: intersection of the sphere
		 * with plans going through (0,0,0) if we are inclusive, the boundaries
		 * are the intersection of the sphere with plains pushed outward by
		 * sin(offset)
		 */
//		double temp = 0.;
		boolean found = false;
		for (iz = irmin; iz <= irmax; iz++) {
			found = false;
			if (iz <= nside - 1) { // North polar cap
				z = 1.0 - iz * iz * dth1;
			} else if (iz <= 3 * nside) { // tropical band + equator
				z = (2.0 * nside - iz) * dth2;
			} else {
				z = -1.0 + (4.0 * nside - iz) * (4.0 * nside - iz) * dth1;
			}

			/* computes the 3 intervals described by the 3 great circles */
			st = Math.sqrt((1.0 - z) * (1.0 + z));
			tgth = z / st; // cotan(theta_ring)
			for (int j = 0; j < 3; j++) {
				dc[j] = tgthi[j] * tgth - sdet * sin_off
						/ ((sto[j] + 1.0e-30) * st);

			}
			for (int k = 0; k < 3; k++) {
				if (dc[k] * sdet <= -1.0) { // the whole iso-latitude ring is on
					// right side of the great circle
					dom[k][0] = 0.0;
					dom[k][1] = TWOPI;
				} else if (dc[k] * sdet >= 1.0) { // all on the wrong side
					dom[k][0] = -1.000001 * (k + 1);
					dom[k][1] = -1.0 * (k + 1);
				} else { // some is good some is bad
					phi_neg = phi0i[k] - (Math.acos(dc[k]) * sdet);
					phi_pos = phi0i[k] + (Math.acos(dc[k]) * sdet);
					//					
					if (phi_pos < 0.)
						phi_pos += TWOPI;
					if (phi_neg < 0.)
						phi_neg += TWOPI;

					//

					dom[k][0] = BitManipulation.MODULO(phi_neg, TWOPI);
					dom[k][1] = BitManipulation.MODULO(phi_pos, TWOPI);

				}
				//

			}
			/* identify the intersections (0,1,2 or 3) of the 3 intervals */

			dom12 = intrs_intrv(dom[0], dom[1]);
			n12 = dom12.length / 2;
			if (n12 != 0) {
				if (n12 == 1) {
					dom123a = intrs_intrv(dom[2], dom12);
					n123a = dom123a.length / 2;

					if (n123a == 0)
						found = true;
					if (!found) {
						for (int l = 0; l < dom123a.length; l++) {
							alldom[l] = dom123a[l];
						}

						ndom = n123a; // 1 or 2
					}
				}
				if (!found) {
					if (n12 == 2) {
						double[] tmp = { dom12[0], dom12[1] };
						dom123a = intrs_intrv(dom[2], tmp);
						double[] tmp1 = { dom12[2], dom12[3] };
						dom123b = intrs_intrv(dom[2], tmp1);
						n123a = dom123a.length / 2;
						n123b = dom123b.length / 2;
						ndom = n123a + n123b; // 0, 1, 2 or 3

						if (ndom == 0)
							found = true;
						if (!found) {
							if (n123a != 0) {
								for (int l = 0; l < 2 * n123a; l++) {
									alldom[l] = dom123a[l];
								}
							}
							if (n123b != 0) {
								for (int l = 0; l < 2 * n123b; l++) {
									alldom[(int) (l + 2 * n123a)] = dom123b[l];
								}
							}
							if (ndom > 3) {
								throw new HealpixException(SID
										+ ": too many intervals found");
							}
						}
					}
				}
				if (!found) {
					for (idom = 0; idom < ndom; idom++) {
						a_i = alldom[(int) (2 * idom)];
						b_i = alldom[(int) (2 * idom + 1)];
						phi0 = (a_i + b_i) / 2.0;
						dphiring = Math.abs(b_i - a_i) / 2.0;

						if (dphiring < 0.0) {
							phi0 += PI;
							dphiring += PI;
						}

						/* finds pixels in the triangle on that ring */
						listir = InRing(nside, iz, phi0, dphiring, do_nest);
						res.addAll(listir);

					}
				}
			}

		}
		return res;
	}

	/**
	 * computes the intersection di of 2 intervals d1 (= [a1,b1])
	 * and d2 (= [a2,b2]) on the periodic domain (=[A,B] where A and B
	 * arbitrary) ni is the resulting number of intervals (0,1, or 2) if a1 <b1
	 * then d1 = {x |a1 <= x <= b1} if a1>b1 then d1 = {x | a1 <=x <= B U A <=x
	 * <=b1}
	 * 
	 * @param d1 double[] first interval
	 * @param d2 double[] second interval
	 * @return double[] one or two intervals intersections
	 */
	public double[] intrs_intrv(double[] d1, double[] d2) {
		double[] res;
		double epsilon = 1.0e-10;
//		double temp = 0.;
//		int ni;
		double[] dk;
		double[] di = { 0. };
		int ik = 0;
		boolean tr12, tr21, tr34, tr43, tr13, tr31, tr24, tr42, tr14, tr32;
		/*                                             */

		tr12 = (d1[0] < d1[1] + epsilon);
		tr21 = !tr12; // d1[0] >= d1[1]
		tr34 = (d2[0] < d2[1] + epsilon);
		tr43 = !tr34; // d2[0]>d2[1]
		tr13 = (d1[0] < d2[0] + epsilon); //  d2[0] can be in interval
		tr31 = !tr13; // d1[0] in longerval
		tr24 = (d1[1] < d2[1] + epsilon); // d1[1] upper limit
		tr42 = !tr24; // d2[1] upper limit
		tr14 = (d1[0] < d2[1] + epsilon); // d1[0] in interval
		tr32 = (d2[0] < d1[1] + epsilon); // d2[0] in interval

		ik = 0;
		dk = new double[] { -1.0e9, -1.0e9, -1.0e9, -1.0e9 };
		/* d1[0] lower limit case 1 */
		if ((tr34 && tr31 && tr14) || (tr43 && (tr31 || tr14))) {
			ik++; // ik = 1;
			dk[ik - 1] = d1[0]; // a1

		}
		/* d2[0] lower limit case 1 */
		if ((tr12 && tr13 && tr32) || (tr21 && (tr13 || tr32))) {
			ik++; // ik = 1
			dk[ik - 1] = d2[0]; // a2

		}
		/* d1[1] upper limit case 2 */
		if ((tr34 && tr32 && tr24) || (tr43 && (tr32 || tr24))) {
			ik++; // ik = 2
			dk[ik - 1] = d1[1]; // b1

		}
		/* d2[1] upper limit case 2 */
		if ((tr12 && tr14 && tr42) || (tr21 && (tr14 || tr42))) {
			ik++; // ik = 2
			dk[ik - 1] = d2[1]; // b2

		}
		di = new double[1];
		di[0] = 0.;
		switch (ik) {

		case 2:
			di = new double[2];

			di[0] = dk[0] - epsilon;
			di[1] = dk[1] + epsilon;
			break;
		case 4:

			di = new double[4];
			di[0] = dk[0] - epsilon;
			di[1] = dk[3] + epsilon;
			di[2] = dk[1] - epsilon;
			di[3] = dk[2] + epsilon;
			break;
		}
		res = di;

		return res;
	}

	/**
	 * an obsolete method. Use query_disc instead.
	 * 
	 * @param nside
	 * @param vector0
	 * @param radius
	 * @return - ArrayList of long
	 */
	public ArrayList getDisc_ring(long nside, Vector3d vector0, double radius) {
		ArrayList res;
		int nest = 0;
		int inclusive = 0;
		res = query_disc(nside, vector0, radius, nest, inclusive);
		return res;
	}

	/**
	 * generates in the RING or NESTED scheme all pixels that lays within an
	 * angular distance Radius of the center.
	 * 
	 * @param nside 
	 *            long map resolution
	 * @param vector 
	 *            Vector3d pointing to the disc center
	 * @param radius 
	 *            double angular radius of the disc (in RADIAN )
	 * @param nest 
	 *            int 0 (default) if output is in RING scheme, if set to 1
	 *            output is in NESTED
	 * @param inclusive 
	 *            int 0 (default) only pixels whose centers lay in the disc
	 *            are listed, if set to 1, all pixels overlapping the disc
	 *            are listed. In the inclusive mode the radius is increased by half the pixel size.
	 *            In this case most probably all neighbor pixels will be listed even with very small
	 *            radius.
	 *            In case of exclusive search and very small radius when the disc lays completely
	 *            inside a pixel the pixel number is returned using vector2pix method.
	 * @return  ArrayList of pixel numbers
	 * 
	 * calls: RingNum(nside, ir) InRing(nside, iz, phi0, dphi,nest) vector2pix(nside,ipix)
	 */
	public ArrayList query_disc(long nside, Vector3d vector, double radius,
			int nest, int inclusive)  {
		ArrayList res = new ArrayList();
		long irmin, irmax, iz, ip, nir, npix, ilist;

//		double norm_vect0;
		double x0, y0, z0, radius_eff;
		double a, b, c, cosang;
		double dth1, dth2;
		double phi0, cosphi0, cosdphi, dphi;
		double rlat0, rlat1, rlat2, zmin, zmax, z;
//		long status, list_size, nlost;
		boolean do_inclusive = false;
		boolean do_nest = false;
		String SID = "QUERY_DISC";
		/*                             */
		long nsidesq = nside * nside;
		npix = 12 * nsidesq;
		double pixres = PixRes(nside);      // in arc seconds
//		double halfPix = 0.5*Math.toRadians(pixres/3600.);  // in radians
		double halfPix = PI / (4.0 * nside);
//		System.out.println("0.5 pixel size ="+halfPix);


		if (radius < 0.0 || radius > PI) {
			throw new IllegalArgumentException(SID 
					+ ": angular radius is in RADIAN and should be in [0,pi]");
		}
		
		if (inclusive == 1)
			do_inclusive = true;
		if (nest == 1)
			do_nest = true;

		dth1 = 1.0 / (3.0 * nside * nside);
		dth2 = 2.0 / (3.0 * nside);

		radius_eff = radius;		
//		if (radius_eff <= halfPix) radius_eff = halfPix;

		/* disc center */
		vector.normalize();
		x0 = vector.x; // norm_vect0;
		y0 = vector.y; // norm_vect0;
		z0 = vector.z; // norm_vect0;
//		System.out.println("x0="+x0+" y0="+y0+" z0="+z0);
		//
		// NK make radius increase a function of altitude
		//
		if (do_inclusive) { radius_eff += (halfPix + Math.abs(z0)*halfPix);}   // increase radius for incluzive search
//		System.out.println("effective radius="+radius_eff);
		//		if (do_inclusive) { radius_eff += halfPix;}  
		cosang = Math.cos(radius_eff);
		phi0 = 0.0;
		dphi = 0.0;
		if (x0 != 0. || y0 != 0.)
			phi0 = BitManipulation.MODULO(Math.atan2(y0, x0) + TWOPI, TWOPI);  // in [0, 2pi]
			cosphi0 = Math.cos(phi0);
//			System.out.println("phi0="+phi0+" cosphi0="+cosphi0);
		a = x0 * x0 + y0 * y0;
		/* coordinate z of highest and lowest points in the disc */
		rlat0 = Math.asin(z0); // latitude in RAD of the center
		rlat1 = rlat0 + radius_eff;
		rlat2 = rlat0 - radius_eff;
//		System.out.println("rlat0="+rlat0+" rlat1="+rlat1+" rlat2="+rlat2);
		//
		if (rlat1 >= HALFPI) {
			zmax = 1.0;
		} else {
			zmax = Math.sin(rlat1);
		}
		irmin = RingNum(nside, zmax);
		irmin = Math.max(1, irmin - 1); // start from a higher point to be safe
		if (rlat2 <= -HALFPI) {
			zmin = -1.0;
		} else {
			zmin = Math.sin(rlat2);
		}
		irmax = RingNum(nside, zmin);
		irmax = Math.min(4 * nside - 1, irmax + 1); // go down to a lower point
//		System.out.println(" irmax="+irmax+" irmin="+irmin);
		ilist = -1;

        

        //
		/* loop on ring number */
		for (iz = irmin; iz <= irmax; iz++) {
			if (iz <= nside - 1) { // north polar cap
				z = 1.0 - iz * iz * dth1;
			} else if (iz <= 3 * nside) { // tropical band + equator
				z = (2.0 * nside - iz) * dth2;
			} else {
				z = -1.0 + (4.0 * nside - iz) * (4.0 * nside - iz) * dth1;
			}
			/* find phi range in the disc for each z */
			b = cosang - z * z0;
			c = 1.0 - z * z;
			cosdphi = b / Math.sqrt(a * c);
			long done = 0;

			if (Math.abs(x0) <= 1.0e-12 && Math.abs(y0) <= 1.0e-12) {
				cosdphi = -1.0;
				dphi = PI;
				done = 1;
			}
			if (done == 0) {
				if (Math.abs(cosdphi) <= 1.0) {
					dphi = Math.acos(cosdphi); // in [0,pi]
				} else {
					if (cosphi0 >= cosdphi) {
						dphi = PI; // all the pixels at this elevation are in
						// the disc
					} else {
						done = 2; // out of the disc
					}
				}

			}
			if (done < 2) { // pixels in disc
				/* find pixels in the disc */
//				System.out.println("iz="+iz+" phi="+phi0+" dphi="+dphi);
				ArrayList listir = InRing(nside, iz, phi0, dphi, do_nest);
//				System.out.println("ir"+iz);

				res.addAll(listir);				
			}

		}
//
// if  radius less than pixel size check that the pixel number is in the list
//  and add one if it is missing.
//

		long pixel = 0;
		if ( pixres > Math.toDegrees(radius)/3600. ) {
//		    System.out.println("case with r < pix. size");
			if (do_nest) {
				pixel = vect2pix_nest(nside,vector);
			} else {
				pixel = vect2pix_ring(nside,vector);
			}
			if (!res.contains(new Long(pixel))) 
			               res.add(new Long(pixel));
		}
		return res;
	}

	/**
	 * generates in the RING or NESTED scheme all pixels that lays at an
	 * angular distance Radius of the center.
	 * 
	 * @param nside 
	 *            long map resolution
	 * @param vector 
	 *            Vector3d pointing to the ring center
	 * @param radius 
	 *            double angular radius of the ring (in RADIAN )
	 * @param nest 
	 *            int 0 (default) if output is in RING scheme, if set to 1
	 *            output is in NESTED
	 * @param inclusive 
	 *            int 0 (default) only pixels whose centers are crossed by the ring
	 *            are listed, if set to 1, all pixels crossed or touched by the ring
	 *            are listed. In the inclusive mode the radius is increased by half the pixel size.
	 *            In this case most probably all neighbor pixels will be listed even with very small
	 *            radius.
	 *            In case of exclusive search and very small radius when the ring lays completely
	 *            inside a pixel the pixel number is returned using vector2pix method.
	 * @return  ArrayList of pixel numbers
	 * 
	 * calls: RingNum(nside, ir) InRing(nside, iz, phi0, dphi,nest) vector2pix(nside,ipix)
	 */
	public ArrayList query_ring(long nside, Vector3d vector, double radius,
			int nest, int inclusive)  {
		ArrayList res = new ArrayList();
		long irmin, irmax, iz, ip, nir, npix, ilist;

//		double norm_vect0;
		double x0, y0, z0, radius_eff;
		double a, b, c, cosang;
		double dth1, dth2;
		double phi0, cosphi0, cosdphi, dphi;
		double rlat0, rlat1, rlat2, zmin, zmax, z;
//		long status, list_size, nlost;
		boolean do_inclusive = false;
		boolean do_nest = false;
		String SID = "QUERY_RING";
		/*                             */
		long nsidesq = nside * nside;
		npix = 12 * nsidesq;
		double pixres = PixRes(nside);      // in arc seconds
//	
		double halfPix = 0.5*Math.sqrt(PI / (3.0*npix)); // in radians
		double dr = 2.0*halfPix;                   // interval in phi
//		System.out.println("0.5 pixel size ="+halfPix);


		if (radius < 0.0 || radius > PI) {
			throw new IllegalArgumentException(SID 
					+ ": angular radius is in RADIAN and should be in [0,pi]");
		}
		
		if (inclusive == 1)
			do_inclusive = true;
		if (nest == 1)
			do_nest = true;

		dth1 = 1.0 / (3.0 * nside * nside);
		dth2 = 2.0 / (3.0 * nside);

		radius_eff = radius;		
//		if (radius_eff <= halfPix) radius_eff = halfPix;

		/* disc center */
		vector.normalize();
		x0 = vector.x; // norm_vect0;
		y0 = vector.y; // norm_vect0;
		z0 = vector.z; // norm_vect0;
//		System.out.println("x0="+x0+" y0="+y0+" z0="+z0);
		//
		// NK make radius increase a function of altitude
		//
		if (do_inclusive) { radius_eff += (halfPix + Math.abs(z0)*halfPix);}   // increase radius for incluzive search
//		System.out.println("effective radius="+radius_eff);
		//		if (do_inclusive) { radius_eff += halfPix;}  
		cosang = Math.cos(radius_eff);
		phi0 = 0.0;
		dphi = 0.0;
		if (x0 != 0. || y0 != 0.)
			phi0 = BitManipulation.MODULO(Math.atan2(y0, x0) + TWOPI, TWOPI);  // in [0, 2pi]
			cosphi0 = Math.cos(phi0);
//			System.out.println("phi0="+phi0+" cosphi0="+cosphi0);
		a = x0 * x0 + y0 * y0;
		/* coordinate z of highest and lowest points in the disc */
		rlat0 = Math.asin(z0); // latitude in RAD of the center
		rlat1 = rlat0 + radius_eff;
		rlat2 = rlat0 - radius_eff;
//		System.out.println("rlat0="+rlat0+" rlat1="+rlat1+" rlat2="+rlat2);
		//
		if (rlat1 >= HALFPI) {
			zmax = 1.0;
		} else {
			zmax = Math.sin(rlat1);
		}
		irmin = RingNum(nside, zmax);
		irmin = Math.max(1, irmin - 1); // start from a higher point to be safe
		if (rlat2 <= -HALFPI) {
			zmin = -1.0;
		} else {
			zmin = Math.sin(rlat2);
		}
		irmax = RingNum(nside, zmin);
		irmax = Math.min(4 * nside - 1, irmax + 1); // go down to a lower point
		System.out.println(" irmax="+irmax+" irmin="+irmin);
		ilist = -1;

        

        //
		/* loop on ring number */
		for (iz = irmin; iz <= irmax; iz++) {
			if (iz <= nside - 1) { // north polar cap
				z = 1.0 - iz * iz * dth1;
			} else if (iz <= 3 * nside) { // tropical band + equator
				z = (2.0 * nside - iz) * dth2;
			} else {
				z = -1.0 + (4.0 * nside - iz) * (4.0 * nside - iz) * dth1;
			}
			/* find phi range in the disc for each z */
			b = cosang - z * z0;
			c = 1.0 - z * z;
			cosdphi = b / Math.sqrt(a * c);
			int done = 0;

			if (Math.abs(x0) <= 1.0e-12 && Math.abs(y0) <= 1.0e-12) {
				cosdphi = -1.0;
				dphi = PI;
				done = 1;
			}
			
			if (done == 0) {
				if (Math.abs(cosdphi) <= 1.0) {
					dphi = Math.acos(cosdphi); // in [0,pi]
				} else {
					if (cosphi0 >= cosdphi) {
						dphi = PI; // 
					} else {
						done = 2; // out of the ring
					}
				}

			}
			System.out.println(" done="+done);
			if (done < 2) { // pixels in ring
				/* find pixels in the ring */
				System.out.println("iz="+iz+" phi="+phi0+" dphi="+dphi);
				double phi01 = phi0 - 0.5*dphi;
				double phi02 = phi0 + 0.5*dphi;
				ArrayList listir = null;

//					System.out.println("iz="+iz +" phi01="+phi01+" phi02="+phi02+" dr="+dr);
					listir = InRing(nside, iz, phi0, dphi, do_nest);
					if (listir.size()<=1) {
					res.addAll(listir);
					} else {
						res.add(listir.get(0));
						res.add(listir.get(listir.size()-1));
					}
//				System.out.println("ir"+iz);
//					listir = InRing(nside, iz, phi02, dr, do_nest);
//					res.addAll(listir);
				
			}

		}
//

		return res;
	}
	
	/**
	 * renders theta and phi coordinates of the nominal pixel center for the
	 * pixel number ipix (RING scheme) given the map resolution parameter nside
	 * 
	 * @param nside 
	 *            long map resolution
	 * @param ipix 
	 *            long pixel number
	 * @return double[] theta,phi
	 */
	public double[] pix2ang_ring(long nside, long ipix)  {
		double[] res = { 0., 0. };
		long nl2, nl4, npix, ncap, iring, iphi, ip, ipix1;
		double fodd, hip, fihip, theta, phi;
		String SID = "pix2ang_ring:";
		/*                            */
		if (nside < 1 || nside > ns_max) {
			throw new IllegalArgumentException(SID + " Nside should be power of 2 >0 and < "+ns_max);
		}
		long nsidesq = nside * nside;
		npix = 12 * nsidesq; // total number of pixels
		if (ipix < 0 || ipix > npix - 1) {
			throw new IllegalArgumentException(SID + " ipix out of range calculated from nside");
		}
		ipix1 = ipix + 1; //  in [1, npix]
		nl2 = 2 * nside;
		nl4 = 4 * nside;
		ncap = 2 * nside * (nside - 1); // points in each polar cap, =0 for
		// nside =1

		if (ipix1 <= ncap) { // North polar cap
			hip = ipix1 / 2.0;
			fihip = (long) hip; // get integer part of hip
			iring = (long) (Math.sqrt(hip - Math.sqrt(fihip))) + 1; // counted from north
			                                                       // pole
			iphi = ipix1 - 2 * iring * (iring - 1);
			theta = Math.acos(1.0 - iring * iring / (3.0 * nsidesq));
			phi = ((double)iphi - 0.5) * PI / (2.0 * iring);


		} else if (ipix1 <= nl2 * (5 * nside + 1)) { // equatorial region
			ip = ipix1 - ncap - 1;
			iring = (long) (ip / nl4) + nside; // counted from North pole
			iphi = (long) BitManipulation.MODULO(ip, nl4) + 1;
			fodd = 0.5 * (1. + BitManipulation.MODULO(iring + nside, 2)); // 1 if iring+nside
			                                                 // is odd, 1/2 otherwise
			theta = Math.acos((nl2 - iring) / (1.5 * nside));
			phi = ((double)iphi - fodd) * PI / (2.0 * nside);

		} else { // South pole cap
			ip = npix - ipix1 + 1;
			hip = ip / 2.0;
			fihip = (long) hip;
			iring = (long) (Math.sqrt(hip - Math.sqrt(fihip))) + 1; // counted from South
			                                                       // pole
			iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1));
			theta = Math.acos(-1.0 + iring * iring / (3.0 * nsidesq));
			phi = ((double)iphi - 0.5) * PI / (2.0 * iring);

		}
		res[0] = theta;
		res[1] = phi;
		return res;
	}

	/**
	 * returns the vector pointing in the center of the pixel ipix. The vector
	 * is calculated by makePix2Vect_ring method
	 * 
	 * @param nside map resolution
	 * @param ipix pixel number
	 * @return Vector3d
	 */
	public Vector3d pix2vect_ring(long nside, long ipix)  {

		PixInfo pixInfo = makePix2Vect_ring(nside, ipix);
		Vector3d res = new Vector3d(pixInfo.pixVect);
		return res;
	}

	/**
	 * returns double [][] with coordinates of the pixel corners. The array is
	 * calculated by makePix2Vect_ring method
	 * 
	 * @param nside map resolution
	 * @param ipix pixel number
	 * @return  double[][] list of vertex coordinates
	 */
	public double[][] pix2vertex_ring(long nside, long ipix)  {
		double[][] res;
		PixInfo pixinfo = makePix2Vect_ring(nside, ipix);
		res = pixinfo.pixVertex;
		return res;
	}

	/**
	 * renders vector (x,y,z) coordinates of the nominal pixel center for pixel
	 * ipix (RING scheme) given the map resolution parameter nside. It also
	 * calculates (x,y,z) positions of the four vertices in order N,W,S,E.
	 * These results are returned in a PixInfo object.
	 * Those can be used using pix2Vect_ring and pix2vert_ring methods
	 * 
	 * @param nside 
	 *            long map resolution
	 * @param ipix 
	 *            pixel number
	 * @return  result object
	 */
	private PixInfo makePix2Vect_ring(long nside, long ipix)  {
		long nl2;
		Vector3d pixVect = new Vector3d(0., 0., 0.);
		double[][] pixVertex = new double[3][4];
                
        long nl4;
        long iring, iphi, ip, ipix1;
        long npix,ncap;
		double phi_nv, phi_wv, phi_sv, phi_ev;
		double z_nv, z_sv, sth_nv, sth_sv, hdelta_phi;
		double fact1, fact2, fodd, hip, fihip, z, sth, phi;
		long iphi_mod;
        long iphi_rat;
//		boolean do_vertex = true;
		long nsidesq = nside * nside;
		String SID = " Pix2Vect_ring:";
		/*                                 */
		if (nside < 1 || nside > ns_max) {
			throw new IllegalArgumentException(SID + " Nside should be power of 2 >0 and < "+ns_max);
		}

		npix = 12 * nsidesq;
		if (ipix < 0 || ipix > npix - 1) {
			throw new IllegalArgumentException(SID + " ipix out of range calculated from nside");
		}

		ipix1 = ipix + 1; //  in [1, npix]
		nl2 = 2 * nside;
		nl4 = 4 * nside;
		ncap = 2 * nside * (nside - 1); // points in each polar cap
		fact1 = 1.5 * nside;
		fact2 = 3.0 * nsidesq;
		phi_nv = 0.0;
		phi_sv = 0.0;
		if (ipix1 <= ncap) { // north polar cap
			hip = ipix1 / 2.0;
			fihip = (long) hip;
			iring = (long) (Math.sqrt(hip - Math.sqrt(fihip))) + 1; // counted from north
			                                                       // pole
			iphi = ipix1 - 2 * iring * (iring - 1);
			z = 1.0 - iring * iring / fact2;
			phi = (iphi - 0.5) * PI / (2.0 * iring);

			hdelta_phi = PI / (4.0 * iring); // half pixel width
			z_nv = 1.0 - (iring - 1) * (iring - 1) / fact2;
			z_sv = 1.0 - (iring + 1) * (iring + 1) / fact2;
			iphi_mod = (long) BitManipulation.MODULO(iphi - 1, iring); // in [0,1,...,iring-1]
			iphi_rat = (iphi - 1) / iring; // in [0,1,2,3]
			if (iring > 1)
				phi_nv = HALFPI * (iphi_rat + iphi_mod / (iring - 1.0));
			phi_sv = HALFPI * (iphi_rat + (iphi_mod + 1.0) / (iring + 1.0));
		} else if (ipix1 <= nl2 * (5 * nside + 1)) { // equatorial region
			ip =  (ipix1 - ncap - 1);
			iring = (long) (ip / nl4) + nside; // counted from North pole
			iphi = (long) BitManipulation.MODULO(ip, nl4) + 1;
			fodd = 0.5 * (1. + BitManipulation.MODULO(iring + nside, 2)); // 1 if iring+nside
			                                                 // is odd or 1/2
			z = (nl2 - iring) / fact1;
			phi = (iphi - fodd) * PI / (2.0 * nside);
			hdelta_phi = PI / (4.0 * nside); // half pixel width
			phi_nv = phi;
			phi_sv = phi;
			z_nv = (nl2 - iring + 1) / fact1;
			z_sv = (nl2 - iring - 1) / fact1;
			if (iring == nside) { // nothern transition
				z_nv = 1.0 - (nside - 1) * (nside - 1) / fact2;
				iphi_mod = (long) BitManipulation.MODULO(iphi - 1, nside); // in [0,1,...,nside-1]
				iphi_rat = (iphi - 1) / nside; // in [0,1,2,3]
				if (nside > 1)
					phi_nv = HALFPI * (iphi_rat + iphi_mod / (nside - 1.));
			} else if (iring == 3 * nside) { // southern transition
				z_sv = -1.0 + (nside - 1) * (nside - 1) / fact2;
				iphi_mod = (long) BitManipulation.MODULO(iphi - 1, nside); // in [0,1,... iring-1]
				iphi_rat = (iphi - 1) / nside; // in [0,1,2,3]
				if (nside > 1)
					phi_sv = HALFPI * (iphi_rat + iphi_mod / (nside - 1.0));
			}

		} else { // South polar cap
			ip = npix - ipix1 + 1;
			hip = ip / 2.0;
			fihip = (long) hip;
			iring = (long) (Math.sqrt(hip - Math.sqrt(fihip))) + 1; // counted from South
			                                                       // pole
			iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1));
			z = -1.0 + iring * iring / fact2;
			phi = (iphi - 0.5) * PI / (2.0 * iring);
			hdelta_phi = PI / (4.0 * iring); // half pixel width
			z_nv = -1.0 + (iring + 1) * (iring + 1) / fact2;
			z_sv = -1.0 + (iring - 1) * (iring - 1) / fact2;
			iphi_mod = (long) BitManipulation.MODULO(iphi - 1, iring); // in [0,1,...,iring-1]
			iphi_rat = (iphi - 1) / iring; // in [0,1,2,3]
			phi_nv = HALFPI * (iphi_rat + (iphi_mod + 1) / (iring + 1.0));
			if (iring > 1)
				phi_sv = HALFPI * (iphi_rat + iphi_mod / (iring - 1.0));

		}
		/* pixel center */
		sth = Math.sqrt((1.0 - z) * (1.0 + z));
		pixVect.x = sth * Math.cos(phi);
		pixVect.y = sth * Math.sin(phi);
		pixVect.z = z;
		pixVect = new Vector3d(sth * Math.cos(phi), sth * Math.sin(phi), z);
		/* west vertex */
		phi_wv = phi - hdelta_phi;
		pixVertex[0][1] = sth * Math.cos(phi_wv);
		pixVertex[1][1] = sth * Math.sin(phi_wv);
		pixVertex[2][1] = z;
		/* east vertex */
		phi_ev = phi + hdelta_phi;
		pixVertex[0][3] = sth * Math.cos(phi_ev);
		pixVertex[1][3] = sth * Math.sin(phi_ev);
		pixVertex[2][3] = z;
		/* north vertex */
		sth_nv = Math.sqrt((1.0 - z_nv) * (1.0 + z_nv));
		pixVertex[0][0] = sth_nv * Math.cos(phi_nv);
		pixVertex[1][0] = sth_nv * Math.sin(phi_nv);
		pixVertex[2][0] = z_nv;
		/* south vertex */
		sth_sv = Math.sqrt((1.0 - z_sv) * (1.0 + z_sv));
		pixVertex[0][2] = sth_sv * Math.cos(phi_sv);
		pixVertex[1][2] = sth_sv * Math.sin(phi_sv);
		pixVertex[2][2] = z_sv;
		return new PixInfo(pixVect, pixVertex);
	}

	/**
	 * renders the pixel number ipix (RING scheme) for a pixel which contains a
	 * point with coordinates theta and phi, given the map resolution parameter
	 * nside.
	 * 
	 * @param nside 
	 *            long map resolution parameter
	 * @param theta 
	 *            double theta
	 * @param phi -
	 *            double phi
	 * @return  long ipix
	 */
	public long ang2pix_ring(long nside, double theta, double phi) {
		long nl4;
        long jp, jm, kshift;
        long ip;
        long ir;
		double z, za, tt, tp, tmp;
		long pix = 0;
		long ipix1;
		long nl2,  ncap, npix;
		String SID = "ang2pix_ring:";
		/*                                       */
		if (nside < 1 || nside > ns_max) {
			throw new IllegalArgumentException(SID + " Nside should be power of 2 >0 and < "+ns_max);
		}
		if (theta < 0.0 || theta > PI) {
			throw new IllegalArgumentException(SID + " Theta out of range [0,pi]");
		}
		
		z = Math.cos(theta);
		za = Math.abs(z);



		if (phi >= TWOPI)  phi = phi -TWOPI ;

		if (phi < 0.)
			phi =phi + TWOPI; //  phi in [0, 2pi]
		tt = phi / HALFPI; // in [0,4]
//		tt = BitManipulation.MODULO(phi, TWOPI) / HALFPI; // in [0,4]
		nl2 = 2 * nside;
		nl4 = 4 * nside;
		ncap = nl2 * (nside - 1); // number of pixels in the north polar cap
		npix = 12 * nside * nside;
		if (za < twothird) { // equatorial region
			jp = (long) (nside * (0.5 + tt - 0.75 * z)); // index of ascending
			// edge line
			jm = (long) (nside * (0.5 + tt + 0.75 * z)); // index of descending
			// edge line

			ir = nside + 1 + jp - jm; // in [1,2n+1]
			kshift = 0;
			if ((long) BitManipulation.MODULO(ir, 2) == 0)
				kshift = 1; // 1 if ir even, 0 otherwise
			ip = (long) ((jp + jm - nside + kshift + 1) / 2) + 1; // in [1,4n]
			if (ip > nl4) ip = ip - nl4;
			ipix1 = ncap + nl4 * (ir - 1) + ip;
			
		} else { // North and South polar caps
			tp = tt - (long) tt;
			tmp = Math.sqrt(3.0 * (1.0 - za));
			jp = (long) (nside * tp * tmp); // increasing edge line index
			jm = (long) (nside * (1.0 - tp) * tmp); // decreasing edge index

			ir = jp + jm + 1; // ring number counted from closest pole
			ip = (long) (tt * ir) + 1; // in [1,4*ir]
			if (ip > 4 * ir)
				ip = ip - 4 * ir;

			ipix1 = 2 * ir * (ir - 1) + ip;
			if (z <= 0.0)
				ipix1 = npix - 2 * ir * (ir + 1) + ip;
						
		}
		pix = ipix1 - 1; // in [0, npix-1]
		

		return pix;
	}

	/**
	 * renders the pixel number ipix (RING scheme) for a pixel which contains a
	 * point on a sphere at coordinate vector (x,y,z), given the map resolution
	 * parameter nside
	 * 
	 * @param nside 
	 *            long map resolution
	 * @param vector 
	 *            Vector3d of the point coordinates
	 * @return  long pixel number
	 * @throws IllegalArgumentException
	 */
	public long vect2pix_ring(long nside, Vector3d vector)  {
		long res = 0;
		long nl2, nl4, ncap, npix, jp, jm, ipix1;
		double z, za, tt, tp, tmp, dnorm, phi;
		long ir, ip, kshift;
		String SID = " vect2pix_ring:";
		/*                                      */
		if (nside < 1 || nside > ns_max) {
			throw new IllegalArgumentException(SID + " Nside should be power of 2 >0 and < "+ns_max);
		}
		dnorm = vector.length();
		z = vector.z / dnorm;
		phi = 0.;
		if (vector.x != 0. || vector.y != 0.)
			phi = Math.atan2(vector.y, vector.x); // phi in [-pi,pi]
		za = Math.abs(z);
		if (phi < 0.)
			phi += TWOPI; //  phi in [0, 2pi]
		tt = phi / HALFPI; // in [0,4]

		nl2 = 2 * nside;
		nl4 = 4 * nside;
		ncap = nl2 * (nside - 1); // number of pixels in the north polar cap
		npix = 12 * nside * nside;
		if (za < twothird) { // equatorial region
			jp = (long) (nside * (0.5 + tt - 0.75 * z)); // index of ascending
			// edge line
			jm = (long) (nside * (0.5 + tt + 0.75 * z)); // index of descending
			// edge line

			ir = nside + 1 + jp - jm; // in [1,2n+1]
			kshift = 0;
			if ((long) BitManipulation.MODULO(ir, 2) == 0)
				kshift = 1; // 1 if ir even, 0 otherwise
			ip = (long) ((jp + jm - nside + kshift + 1) / 2) + 1; // in [1,4n]
			ipix1 = ncap + nl4 * (ir - 1) + ip;
		} else { // North and South polar caps
			tp = tt - (long) tt;
			tmp = Math.sqrt(3.0 * (1.0 - za));
			jp = (long) (nside * tp * tmp); // increasing edge line index
			jm = (long) (nside * (1.0 - tp) * tmp); // decreasing edge index

			ir = jp + jm + 1; // ring number counted from closest pole
			ip = (long) (tt * ir) + 1; // in [1,4*ir]
			if (ip > 4 * ir)
				ip = ip - 4 * ir;

			ipix1 = 2 * ir * (ir - 1) + ip;
			if (z <= 0.0)
				ipix1 = npix - 2 * ir * (ir + 1) + ip;
		}
		res = ipix1 - 1; // in [0, npix-1]
		return res;
	}

	/**
	 * 
	 * Renders theta and phi coordinates of the normal pixel center for the
	 * pixel number ipix (NESTED scheme) given the map resolution parameter
	 * nside.
	 * 
	 * @param nside 
	 *            map resolution parameter - long
	 * @param ipix 
	 *            long pixel number 
	 * @return double[] (theta, phi)
	 * @throws IllegalArgumentException
	 */
	public double[] pix2ang_nest(long nside, long ipix)  {
		double[] res = new double[2];
		double theta = 0.;
		double phi = 0.;
		long npix, npface, ipf, ip_low, ip_trunc, ip_med, ip_hi;
		long jrt, jr, nr, jpt, jp, kshift, nl4, ix, iy, face_num;
		double z, fn, fact1, fact2;
		String SID = "pix2ang_nest:";
		// coordinate of the lowest corner of each face
		long[] jrll = { 0, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 }; // in units of
		// nside
		long[] jpll = { 0, 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 }; // in units of
		// nside/2
		/*                                                              */
		if (nside < 1 || nside > ns_max) {
			throw new IllegalArgumentException(SID + " Nside should be power of 2 >0 and < "+ns_max);
		}
		long nsidesq = nside * nside;
		npix = 12 * nsidesq;
		if (ipix < 0 || ipix > npix - 1) {
			throw new IllegalArgumentException(SID + " ipix out of range calculated from nside");
		}
		if (pix2x[xmax-1] <= 0)
			mk_pix2xy();
		fn = 1.*nside;
		fact1 = 1.0 / (3.0 * fn * fn);
		fact2 = 2.0 / (3.0 * fn);
		nl4 = 4 * nside;
		/* findes the face, and the number in the face */
		npface = nside * nside;
		face_num = ipix / npface; // face number [0,11]
		ipf = (long) BitManipulation.MODULO(ipix, npface); // pixel in the face [0, npface-1]
		/*
		 * finds x,y on the face (starting from the lowest corner) from pixel
		 * number
		 */
		ip_low = (long) BitManipulation.MODULO(ipf, pixmax);   // content of the last 18 bits
		ip_trunc = ipf / pixmax;                  // trancation of the last 18 bits
		ip_med = (long) BitManipulation.MODULO(ip_trunc, pixmax); // content of the next 18 bits
		ip_hi = ip_trunc / pixmax;                // content of the high wait 18 bits

		ix = pixmax * pix2x[(int)ip_hi] + xmid * pix2x[(int)ip_med] + pix2x[(int) ip_low];
		iy = pixmax * pix2y[(int)ip_hi] + xmid * pix2y[(int)ip_med] + pix2y[(int)ip_low];
		/* transform these in (horizontal, vertical) coordinates */
		jrt = ix + iy; // [0,2*(nside-1)]
		jpt = ix - iy; // [ -nside+1, nside -1]
		/* computes the z coordinate on the sphere */
		jr = jrll[(int) (face_num + 1)] * nside - jrt - 1; // ring number in [1,
		// 4*nside-1]

		nr = nside; // equatorial region (the most frequent )
		z = (2 * nside - jr) * fact2;
		kshift = (long) BitManipulation.MODULO(jr - nside, 2);
		if (jr < nside) { // north pole region
			nr = jr;
			z = 1.0 - nr * nr * fact1;
			kshift = 0;
		} else if (jr > 3 * nside) { // south pole region
			nr = nl4 - jr;
			z = -1.0 + nr * nr * fact1;
			kshift = 0;
		}
		theta = Math.acos(z);
		/* computes phi coordinate on the sphere, in [0,2pi] */
		jp = (jpll[(int) (face_num + 1)] * nr + jpt + 1 + kshift) / 2;
		if (jp > nl4)
			jp = jp - nl4;
		if (jp < 1)
			jp = jp + nl4;

		phi = (jp - (kshift + 1) * 0.5) * (HALFPI / nr);
		res[0] = theta;
		res[1] = phi;
		return res;
	}

	/**
	 * renders vector (x,y,z) coordinates of the nominal pixel center for the
	 * pixel ipix (NESTED scheme ) given the map resolution parameter nside.
	 * Also calculates the (x,y,z) positions of 4 pixel vertices (corners) in
	 * the order N,W,S,E. These can be get using method pix2vertex_nest.
	 * 
	 * @param nside the map resolution 
	 * @param ipix long pixel number
	 * @return Vector3d
	 * @throws IllegalArgumentException
	 */
	public Vector3d pix2vect_nest(long nside, long ipix)  {

		PixInfo pixinfo = makePix2Vect_Nest(nside, ipix);
		return pixinfo.pixVect;
	}

	/**
	 * renders vector (x,y,z) coordinates of the nominal pixel center for the
	 * pixel ipix (NESTED scheme ) given the map resolution parameter nside.
	 * Also calculates the (x,y,z) positions of 4 pixel vertices (corners) in
	 * the order N,W,S,E.
	 * 
	 * @param nside the map resolution
	 * @param ipix long pixel number
	 * @return double[3][4] 4 sets of vector components
	 * @throws IllegalArgumentException
	 */
	public double[][] pix2vertex_nest(long nside, long ipix) {
		double[][] res = new double[3][4];
		PixInfo pixinfo = makePix2Vect_Nest(nside, ipix);
		double[][] pixVertex = pixinfo.pixVertex;
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 4; j++) {
				res[i][j] = pixVertex[i][j];
			}
		}
		return res;
	}

	/**
	 * renders vector (x,y,z) coordinates of the nominal pixel center for the
	 * pixel ipix (NESTED scheme ) given the map resolution parameter nside. Tis
	 * can be get using method pix2vect_nest Also calculates the (x,y,z)
	 * positions of 4 pixel vertices (corners) in the order N,W,S,E.
	 * The result can be used using method pix2vertex_nest.
	 * @param nside long the map resolution
	 * @param ipix long pixel number
	 * @return result
	 */
	private PixInfo makePix2Vect_Nest(long nside, long ipix)  {
		Vector3d pixVect = new Vector3d(0., 0., 0.);
		double[][] pixVertex = new double[3][4];
		long npix, npface, ipf, ip_low, ip_trunc, ip_med, ip_hi;
		long jrt, jr, nr, jpt, jp, kshift, nl4;
		double z, fn, fact1, fact2, sth, phi;
		long ix, iy, face_num;
		int[] jrll = { 0, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 };
		int[] jpll = { 0, 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 };
		double phi_nv, phi_wv, phi_sv, phi_ev, phi_up, phi_dn;
		double z_nv, z_sv, sth_nv, sth_sv;
		double hdelta_phi;
		long iphi_mod, iphi_rat;
//		boolean do_vertex = true;
		String SID = "Pix2Vect_Nest:";
		z_nv = 0.;
		z_sv = 0.;
		/*                                 */
		if (nside < 1 || nside > ns_max) {
			throw new IllegalArgumentException(SID + " Nside should be power of 2 >0 and < "+ns_max);
		}
		long nsidesq = nside * nside;
		npix = 12 * nsidesq;
		if (ipix < 0 || ipix > npix - 1) {
			throw new IllegalArgumentException(SID + " ipix out of range calculated from nside");
		}
		/* initiates the array for the pixel number -> (x,y) mapping */
		if (pix2x[pixmax-1] <= 0)
			mk_pix2xy();
		fn = nside;
		fact1 = 1.0 / (3.0 * fn * fn);
		fact2 = 2.0 / (3.0 * fn);
		nl4 = 4 * nside;

		/* finds the face and the number in the face */
		npface = nsidesq;
		//		System.out.println("ipix="+ipix+" npface="+npface);
		face_num = ipix / npface; // face number in [0, 11]
		ipf = (long) BitManipulation.MODULO(ipix, npface); // pixel number in the face [0,npface-1]
		/*
		 * finds the x,y on the face (starting from the lowest corner) from the
		 * pixel number
		 */
		ip_low = (long) BitManipulation.MODULO(ipf, pixmax); // last 18 bits
		ip_trunc = ipf / pixmax; // trancateion of the last 18 bits
		ip_med = (long) BitManipulation.MODULO(ip_trunc, pixmax); // next 18 bits
		ip_hi = ip_trunc / pixmax; // high 18 bits

		ix = pixmax * pix2x[(int)ip_hi] + xmid * pix2x[(int)ip_med] + pix2x[(int)ip_low];
		iy = pixmax * pix2y[(int)ip_hi] + xmid * pix2y[(int)ip_med] + pix2y[(int)ip_low];

		/* transform this in (vertical, horizontal) coordinates */
		jrt = ix + iy; //  vertical in [0,2*(nside-1)]
		jpt = ix - iy; //  horizontal in [ -nside+1, nside-1]
		/* computes the z coordinate on the sphere */
		jr = jrll[(int) (face_num + 1)] * nside - jrt - 1; // ring number in
		// [1,4*nside-1]

		nr = nside; // equatorial region (the most frequent )
		z = (2.0 * nside - jr) * fact2;

		kshift = (long) BitManipulation.MODULO(jr - nside, 2);
		z_nv = (2.0 * nside - jr + 1.0) * fact2;
		z_sv = (2.0 * nside - jr - 1.0) * fact2;
		if (jr == nside) { // northen transition
			z_nv = 1.0 - (nside - 1.0) * (nside - 1.0) * fact1;
		} else if (jr == 3 * nside) { // southern transition
			z_sv = -1.0 + (nside - 1.0) * (nside - 1.0) * fact1;
		}

		if (jr < nside) { // north pole region
			nr = jr;
			z = 1.0 - nr * nr * fact1;
			kshift = 0;

			z_nv = 1.0 - (nr - 1) * (nr - 1) * fact1;
			z_sv = 1.0 - (nr + 1) * (nr + 1) * fact1;

		} else if (jr > 3 * nside) { // south pole region
			nr = nl4 - jr;
			z = -1.0 + nr * nr * fact1;
			kshift = 0;

			z_nv = -1.0 + (nr + 1) * (nr + 1) * fact1;
			z_sv = -1.0 + (nr - 1) * (nr - 1) * fact1;

		}
		/* computes the phi coordinate on the sphere, in [0,2pi] */
		jp = (jpll[(int) (face_num + 1)] * nr + jpt + 1 + kshift) / 2; // phi in the ring in
		                                                       // [1,4*nr]
		if (jp > nl4)
			jp = jp - nl4;
		if (jp < 1)
			jp = jp + nl4;

		phi = (jp - (kshift + 1) / 2.) * (HALFPI / nr);

		sth = Math.sqrt((1.0 - z) * (1.0 + z));
		pixVect.x = sth * Math.cos(phi);
		pixVect.y = sth * Math.sin(phi);
		pixVect.z = z;

		phi_nv = phi;
		phi_sv = phi;

		phi_up = 0.0;
		iphi_mod = (long) BitManipulation.MODULO(jp - 1, nr); // in [0,1,...,nr-1]
		iphi_rat = (jp - 1) / nr; // in [0,1,2,3]
		if (nr > 1)
			phi_up = HALFPI * (iphi_rat + iphi_mod / (nr - 1.));
		phi_dn = HALFPI * (iphi_rat + (iphi_mod + 1) / (nr + 1.));

		if (jr < nside) { // north polar cap
			phi_nv = phi_up;
			phi_sv = phi_dn;
		} else if (jr > 3 * nside) { // south polar cap
			phi_nv = phi_dn;
			phi_sv = phi_up;
		} else if (jr == nside) { // north transition
			phi_nv = phi_up;
		} else if (jr == 3 * nside) { // south transition
			phi_sv = phi_up;
		}
		hdelta_phi = PI / (4.0 * nr);
		/* west vertex */
//		phi_wv = phi = hdelta_phi;
		phi_wv = phi - hdelta_phi;
		pixVertex[0][1] = sth * Math.cos(phi_wv);
		pixVertex[1][1] = sth * Math.sin(phi_wv);
		pixVertex[2][1] = z;
		/* east vertex */
		phi_ev = phi + hdelta_phi;
		pixVertex[0][3] = sth * Math.cos(phi_ev);
		pixVertex[1][3] = sth * Math.sin(phi_ev);
		pixVertex[2][3] = z;
		/* north vertex */
		sth_nv = Math.sqrt((1.0 - z_nv) * (1.0 + z_nv));
		pixVertex[0][0] = sth_nv * Math.cos(phi_nv);
		pixVertex[1][0] = sth_nv * Math.sin(phi_nv);
		pixVertex[2][0] = z_nv;
		/* south vertex */
		sth_sv = Math.sqrt((1.0 - z_sv) * (1.0 + z_sv));
		pixVertex[0][2] = sth_sv * Math.cos(phi_sv);
		pixVertex[1][2] = sth_sv * Math.sin(phi_sv);
		pixVertex[2][2] = z_sv;
		return new PixInfo(pixVect, pixVertex);
	}

	/**
	 * renders the pixel number pix (NESTED scheme) for a pixel which contains a
	 * point on a sphere at coordinates theta and phi, given map resolution
	 * parameter nside.
	 * 
	 * The computation is made to the highest resolution available and then
	 * degraded to required resolution by integer division. It makes sure that
	 * the treatment of round-off will be consistent for every resolution.
	 * 
	 * @param nside the map resolution parameter
	 * @param theta double theta coordinate in radians
	 * @param phi double phi coordinate in radians
	 * @return pixel number long 
	 * @throws IllegalArgumentException
	 */
	public long ang2pix_nest(long nside, double theta, double phi) {
		long pix = 0;
//		long ipix1;
		double z, za, tt, tp, tmp;
		long jp, jm, ifp, ifm, face_num, ix, iy, ix_low, ix_hi;
		long iy_low, iy_hi, ipf, ntt;
//		long nl2, nl4, ncap, npix, kshift, ir, ip;
		String SID = "ang2pix_nest:";
		/*                              */
		if (nside < 1 || nside > ns_max) {
			throw new IllegalArgumentException(SID + " Nside should be power of 2 >0 and < "+ns_max);
		}
		if (theta < 0.0 || theta > PI) {
			throw new IllegalArgumentException(SID + " theta is out of range [0.,PI]");
		}
		if (x2pix[xmax-1] <= 0)			
			 mk_xy2pix();
			
		z = Math.cos(theta);
		za = Math.abs(z);


		if (phi < 0.)
			phi += TWOPI; // phi in [0,2pi]
		if (phi >= TWOPI ) phi -= TWOPI;
		tt = BitManipulation.MODULO(phi, TWOPI) / HALFPI; // in [0,4]
//		tt = 2. * phi / PI; // in [0,4]
		if (za <= twothird) { // Equatorial region
			/*
			 * the index of edge lines increases when the longitude = phi goes
			 * up
			 */
			jp = (long) (ns_max * (0.5 + tt - z * 0.75)); // ascending edge line
			// index
			jm = (long) (ns_max * (0.5 + tt + z * 0.75)); // descending edge line
			// index
			/* find the face */
			ifp = jp / ns_max; // in [0,4]
			ifm = jm / ns_max;
			if (ifp == ifm) { // faces 4 to 7
				face_num = (long) BitManipulation.MODULO(ifp, 4) + 4;
			} else if (ifp < ifm) { // (half-) faces 0 to 3
				face_num = (long) BitManipulation.MODULO(ifp, 4);
			} else { // (half-) faces 8 to 11
				face_num = (long) BitManipulation.MODULO(ifm, 4) + 8;
			}

			ix = (long) BitManipulation.MODULO(jm, ns_max);
			iy = ns_max - (long) BitManipulation.MODULO(jp, ns_max) - 1;
		} else { // polar region, za > 2/3
			ntt = (long) tt;
			if (ntt >= 4)
				ntt = 3;
			tp = tt - ntt;
			tmp = Math.sqrt(3.0 * (1.0 - za)); // in [0,1]
			/*
			 * the index of edge lines increases when distance from the closest
			 * pole goes up
			 */
			jp = (long) (ns_max * tp * tmp); // line going toward the pole has
			// phi increases
			jm = (long) (ns_max * (1.0 - tp) * tmp); // that one goes away of the
			// closest pole
			jp = (long) Math.min(ns_max - 1, jp); // for pointss too close to the
			// boundary
			jm = (long) Math.min(ns_max - 1, jm);
			/* finds the face and pixel's (x,y) */
			if (z >= 0) {
				face_num = ntt; // in [0,3]
				ix = ns_max - jm - 1;
				iy = ns_max - jp - 1;
			} else {
				face_num = ntt + 8; // in [8,11]
				ix = jp;
				iy = jm;
			}
		}
		ix_low = (long) BitManipulation.MODULO(ix, xmax);
		ix_hi = ix / xmax;
		iy_low = (long) BitManipulation.MODULO(iy, xmax);
		iy_hi = iy / xmax;

		ipf = (x2pix[(int) (ix_hi + 1)] + y2pix[(int) (iy_hi + 1)]) * (xmax * xmax)
				+ (x2pix[(int) (ix_low + 1)] + y2pix[(int) (iy_low + 1)]);
		ipf = ipf / ((ns_max / nside) * (ns_max / nside)); // in [0, nside**2
		// -1]
		pix = ipf + face_num * nside * nside; // in [0, 12*nside**2 -1]
		return pix;
	}

	/**
	 * make the conversion NEST to RING
	 * 
	 * @param nside the map resolution parameter
	 * @param map Object[] the map in NESTED scheme
	 * @return - Object[] a map in RING scheme
	 * @throws IllegalArgumentException
	 */
	public Object[] convert_nest2ring(long nside, Object[] map)  {
		Object[] res;
		long npix, ipn;
        int ipr;
		npix = 12 * nside * nside;
		res = new Object[(int) npix];
		for (ipn = 0; ipn < npix; ipn++) {
			ipr = (int) nest2ring(nside, ipn);
			res[ipr] = map[(int) ipn];
		}
		return res;
	}

	/**
	 * makes the conversion RING to NEST
	 * 
	 * @param nside 
	 *            long resolution
	 * @param map 
	 *            map in RING
	 * @return  map in NEST
	 * @throws IllegalArgumentException
	 */
	public Object[] convert_ring2nest(long nside, Object[] map)  {
		Object[] res;
		long npix, ipn, ipr;
		npix = 12 * nside * nside;
		res = new Object[(int) npix];
		for (ipr = 0; ipr < npix; ipr++) {
			ipn = ring2nest(nside, ipr);
			res[(int) ipn] = map[(int)ipr];
		}
		return res;
	}

	/**
	 * converts a 8 byte Object map from RING to NESTED and vice versa in place,
	 * ie without allocation a temporary map (Has no reason for Java). This
	 * method is more general but slower than convert_nest2ring.
	 * 
	 * This method is a wrapper for functions ring2nest and nest2ring. Their
	 * names are supplied in the subcall argument.
	 * 
	 * @param subcall 
	 *            String name of the method to use.
	 * @param map 
	 *            Object[] map
	 * @return  resulting Object[] map.
	 * @throws IllegalArgumentException
	 */
	public Object[] convert_inplace_long(String subcall, Object[] map) {
		Object[] res;
		long npix, nside;
		boolean[] check;
		long ilast, i1, i2;
		String SID = "convert_in_place:";
		Object pixbuf1, pixbuf2;
		npix = map.length;
		nside = (long) Math.sqrt(npix / 12.);
		if (nside > ns_max) {
			throw new IllegalArgumentException(SID + " Map is too big");
		}
		check = new boolean[(int) npix];
		for (int i = 0; i < npix; i++)
			check[i] = false;
		ilast = 0; // start from first pixel
		for (int i = 0; i < npix; i++) {
			pixbuf2 = map[(int) ilast];
			i1 = ilast;
			if (subcall.equalsIgnoreCase("ring2nest")) {
				i2 = ring2nest(nside, i1);
			} else {
				i2 = nest2ring(nside, i1);
			}
			while (!check[(int) i2]) {
				pixbuf1 = map[(int) i2];
				map[(int) i2] = pixbuf2;
				pixbuf2 = pixbuf1;
				i1 = i2;
				if (subcall.equalsIgnoreCase("ring2nest")) {
					i2 = ring2nest(nside, i1);
				} else {
					i2 = nest2ring(nside, i1);
				}
			}
			while (!(check[(int) ilast] && (ilast < npix - 1))) {
				ilast++;
			}
		}
		res = map;
		return res;
	}

	/**
	 * returns 7 or 8 neighbours of any pixel in the nested scheme The neighbours
	 * are ordered in the following way: First pixel is the one to the south (
	 * the one west of the south direction is taken for pixels that don't have a
	 * southern neighbour). From then on the neighbors are ordered in the
	 * clockwise direction.
	 * 
	 * @param nside the map resolution
	 * @param ipix long pixel number
	 * @return ArrayList
	 * @throws IllegalArgumentException
	 */
	public ArrayList neighbours_nest(long nside, long ipix)  {
		ArrayList res = new ArrayList();
		long npix, ipf, ipo, ix, ixm, ixp, iy, iym, iyp, ixo, iyo;
		long face_num, other_face;
		long ia, ib, ibp, ibm, ib2,  nsidesq;
        int icase;
		long local_magic1, local_magic2;
		long arb_const = 0;
		long[] ixiy = new long[2];
		long[] ixoiyo = new long[2];
		String SID = "neighbours_nest:";
		/* fill the pixel list with 0 */
		res.add(0, new Long(0));
		res.add(1, new Long(0));
		res.add(2, new Long(0));
		res.add(3, new Long(0));
		res.add(4, new Long(0));
		res.add(5, new Long(0));
		res.add(6, new Long(0));
		res.add(7, new Long(0));
		icase = 0;
		/*                                 */
		if ((nside < 1) || (nside > ns_max)) {
			throw new IllegalArgumentException(SID + " Nside should be power of 2 >0 and < "+ns_max);
		}
		nsidesq = nside * nside;
		npix = 12 * nsidesq; // total number of pixels
		if ((ipix < 0) || (ipix > npix - 1)) {
			throw new IllegalArgumentException(SID + " ipix out of range ");
		}
		if (x2pix[xmax-1] <= 0)
			mk_xy2pix();

		local_magic1 = (nsidesq - 1) / 3;
		local_magic2 = 2 * local_magic1;
		face_num = ipix / nsidesq;
		ipf = (long) BitManipulation.MODULO(ipix, nsidesq); // Pixel number in face
		ixiy = pix2xy_nest(nside, ipf);
		ix = ixiy[0];
		iy = ixiy[1];
		//
		ixm = ixiy[0] - 1;
		ixp = ixiy[0] + 1;
		iym = ixiy[1] - 1;
		iyp = ixiy[1] + 1;

		icase = 0; // inside the face

		/* exclude corners */
		if (ipf == local_magic2 && icase == 0)
			icase = 5; // West corner
		if (ipf == (nsidesq - 1) && icase == 0)
			icase = 6; // North corner
		if (ipf == 0 && icase == 0)
			icase = 7; // South corner
		if (ipf == local_magic1 && icase == 0)
			icase = 8; // East corner

		/* detect edges */
		if ((ipf & local_magic1) == local_magic1 && icase == 0)
			icase = 1; // NorthEast
		if ((ipf & local_magic1) == 0 && icase == 0)
			icase = 2; // SouthWest
		if ((ipf & local_magic2) == local_magic2 && icase == 0)
			icase = 3; // NorthWest
		if ((ipf & local_magic2) == 0 && icase == 0)
			icase = 4; // SouthEast

		/* iside a face */
		if (icase == 0) {
			res.add(0, new Long( xy2pix_nest(nside, ixm, iym, face_num)));
			res.add(1, new Long( xy2pix_nest(nside, ixm, iy, face_num)));
			res.add(2, new Long( xy2pix_nest(nside, ixm, iyp, face_num)));
			res.add(3, new Long( xy2pix_nest(nside, ix, iyp, face_num)));
			res.add(4, new Long( xy2pix_nest(nside, ixp, iyp, face_num)));
			res.add(5, new Long( xy2pix_nest(nside, ixp, iy, face_num)));
			res.add(6, new Long( xy2pix_nest(nside, ixp, iym, face_num)));
			res.add(7, new Long( xy2pix_nest(nside, ix, iym, face_num)));
		}
		/*                 */
		ia = face_num / 4; // in [0,2]
		ib = (long) BitManipulation.MODULO(face_num, 4); // in [0,3]
		ibp = (long) BitManipulation.MODULO(ib + 1, 4);
		ibm = (long) BitManipulation.MODULO(ib + 4 - 1, 4);
		ib2 = (long) BitManipulation.MODULO(ib + 2, 4);

		if (ia == 0) { // North pole region
			switch (icase) {
			case 1: // north-east edge
				other_face = 0 + ibp;
				res.set(0, new Long( xy2pix_nest(nside, ixm, iym, face_num)));
				res.set(1, new Long( xy2pix_nest(nside, ixm, iy, face_num)));
				res.set(2, new Long( xy2pix_nest(nside, ixm, iyp, face_num)));
				res.set(3, new Long( xy2pix_nest(nside, ix, iyp, face_num)));
				res.set(7, new Long( xy2pix_nest(nside, ix, iym, face_num)));
				ipo = (long) BitManipulation.MODULO(BitManipulation.swapLSBMSB( ipf), nsidesq);
				ixoiyo = pix2xy_nest(nside, ipo);
				ixo = ixoiyo[0];
				iyo = ixoiyo[1];
				res.set(4, new Long( xy2pix_nest(nside, ixo + 1, iyo,
						other_face)));
				res.set(5, new Long( (other_face * nsidesq + ipo)));
				res.set(6, new Long( xy2pix_nest(nside, ixo - 1, iyo,
						other_face)));
				break;
			case 2: // SouthWest edge
				other_face = 4 + ib;
				ipo = (long) BitManipulation.MODULO(BitManipulation.invLSB( ipf), nsidesq); // SW-NE flip
				ixoiyo = pix2xy_nest(nside, ipo);
				ixo = ixoiyo[0];
				iyo = ixoiyo[1];
				res.set(0, new Long( xy2pix_nest(nside, ixo, iyo - 1,
						other_face)));
				res.set(1, new Long( (other_face * nsidesq + ipo)));
				res.set(2, new Long( xy2pix_nest(nside, ixo, iyo + 1,
						other_face)));
				res.set(3, new Long( xy2pix_nest(nside, ix, iyp, face_num)));
				res.set(4, new Long( xy2pix_nest(nside, ixp, iyp, face_num)));
				res.set(5, new Long( xy2pix_nest(nside, ixp, iy, face_num)));
				res.set(6, new Long( xy2pix_nest(nside, ixp, iym, face_num)));
				res.set(7, new Long( xy2pix_nest(nside, ix, iym, face_num)));
				break;
			case 3: // NorthWest edge
				other_face = 0 + ibm;
				ipo = (long) BitManipulation.MODULO(BitManipulation.swapLSBMSB( ipf), nsidesq); // E-W flip
				ixoiyo = pix2xy_nest(nside, ipo);
				ixo = ixoiyo[0];
				iyo = ixoiyo[1];
				res.set(0, new Long( xy2pix_nest(nside, ixm, iym, face_num)));
				res.set(1, new Long( xy2pix_nest(nside, ixm, iy, face_num)));
				res.set(2, new Long( xy2pix_nest(nside, ixo, iyo - 1,
						other_face)));
				res.set(3, new Long( (other_face * nsidesq + ipo)));
				res.set(4, new Long( xy2pix_nest(nside, ixo, iyo + 1,
						other_face)));
				res.set(5, new Long( xy2pix_nest(nside, ixp, iy, face_num)));
				res.set(6, new Long( xy2pix_nest(nside, ixp, iym, face_num)));
				res.set(7, new Long( xy2pix_nest(nside, ix, iym, face_num)));
				break;
			case 4: // SouthEast edge
				other_face = 4 + ibp;
				ipo = (long) BitManipulation.MODULO(BitManipulation.invMSB( ipf), nsidesq); // SE-NW flip
				ixoiyo = pix2xy_nest(nside, ipo);
				ixo = ixoiyo[0];
				iyo = ixoiyo[1];
				res.set(0, new Long( xy2pix_nest(nside, ixo - 1, iyo,
						other_face)));
				res.set(1, new Long( xy2pix_nest(nside, ixm, iy, face_num)));
				res.set(2, new Long( xy2pix_nest(nside, ixm, iyp, face_num)));
				res.set(3, new Long( xy2pix_nest(nside, ix, iyp, face_num)));
				res.set(4, new Long( xy2pix_nest(nside, ixp, iyp, face_num)));
				res.set(5, new Long( xy2pix_nest(nside, ixp, iy, face_num)));
				res.set(6, new Long( xy2pix_nest(nside, ixo + 1, iyo,
						other_face)));
				res.set(7, new Long( (other_face * nsidesq + ipo)));
				break;
			case 5: // West corner
				other_face = 4 + ib;
				arb_const = other_face * nsidesq + nsidesq - 1;
				res.set(0, new Long( (arb_const - 2)));
				res.set(1, new Long( arb_const));
				other_face = 0 + ibm;
				arb_const = other_face * nsidesq + local_magic1;
				res.set(2, new Long( arb_const));
				res.set(3, new Long( (arb_const + 2)));
				res.set(4, new Long( (ipix + 1)));
				res.set(5, new Long( (ipix - 1)));
				res.set(6, new Long( (ipix - 2)));
				res.remove(7);
				break;
			case 6: //  North corner
				other_face = 0 + ibm;
				res.set(0, new Long( (ipix - 3)));
				res.set(1, new Long( (ipix - 1)));
				arb_const = other_face * nsidesq + nsidesq - 1;
				res.set(2, new Long( (arb_const - 2)));
				res.set(3, new Long( arb_const));
				other_face = 0 + ib2;
				res.set(4, new Long( (other_face * nsidesq + nsidesq - 1)));
				other_face = 0 + ibp;
				arb_const = other_face * nsidesq + nsidesq - 1;
				res.set(5, new Long( arb_const));
				res.set(6, new Long( (arb_const - 1)));
				res.set(7, new Long( (ipix - 2)));
				break;
			case 7: // South corner
				other_face = 8 + ib;
				res.set(0, new Long( (other_face * nsidesq + nsidesq - 1)));
				other_face = 4 + ib;
				arb_const = other_face * nsidesq + local_magic1;
				res.set(1, new Long( arb_const));
				res.set(2, new Long( (arb_const + 2)));
				res.set(3, new Long( (ipix + 2)));
				res.set(4, new Long( (ipix + 3)));
				res.set(5, new Long( (ipix + 1)));
				other_face = 4 + ibp;
				arb_const = other_face * nsidesq + local_magic2;
				res.set(6, new Long( (arb_const + 1)));
				res.set(7, new Long( arb_const));
				break;
			case 8: // East corner
				other_face = 0 + ibp;
				res.set(1, new Long( (ipix - 1)));
				res.set(2, new Long( (ipix + 1)));
				res.set(3, new Long( (ipix + 2)));
				arb_const = other_face * nsidesq + local_magic2;
				res.set(4, new Long( (arb_const + 1)));
				res.set(5, new Long(( arb_const)));
				other_face = 4 + ibp;
				arb_const = other_face * nsidesq + nsidesq - 1;
				res.set(0, new Long( (arb_const - 1)));
				res.set(6, new Long( arb_const));
				res.remove(7);
				break;
			}
		} else if (ia == 1) { // Equatorial region
			switch (icase) {
			case 1: // north-east edge
				other_face = 0 + ibp;
				res.set(0, new Long( xy2pix_nest(nside, ixm, iym, face_num)));
				res.set(1, new Long( xy2pix_nest(nside, ixm, iy, face_num)));
				res.set(2, new Long( xy2pix_nest(nside, ixm, iyp, face_num)));
				res.set(3, new Long( xy2pix_nest(nside, ix, iyp, face_num)));
				res.set(7, new Long( xy2pix_nest(nside, ix, iym, face_num)));
				ipo = (long) BitManipulation.MODULO(BitManipulation.invLSB( ipf), nsidesq); // NE-SW flip
				ixoiyo = pix2xy_nest(nside, ipo);
				ixo = ixoiyo[0];
				iyo = ixoiyo[1];
				res.set(4, new Long( xy2pix_nest(nside, ixo, iyo + 1,
						other_face)));
				res.set(5, new Long( (other_face * nsidesq + ipo)));
				res.set(6, new Long( xy2pix_nest(nside, ixo, iyo - 1,
						other_face)));
				break;
			case 2: // SouthWest edge
				other_face = 8 + ibm;
				ipo = (long) BitManipulation.MODULO(BitManipulation.invLSB( ipf), nsidesq); // SW-NE flip
				ixoiyo = pix2xy_nest(nside, ipo);
				ixo = ixoiyo[0];
				iyo = ixoiyo[1];
				res.set(0, new Long( xy2pix_nest(nside, ixo, iyo - 1,
						other_face)));
				res.set(1, new Long((other_face * nsidesq + ipo)));
				res.set(2, new Long( xy2pix_nest(nside, ixo, iyo + 1,
						other_face)));
				res.set(3, new Long( xy2pix_nest(nside, ix, iyp, face_num)));
				res.set(4, new Long( xy2pix_nest(nside, ixp, iyp, face_num)));
				res.set(5, new Long( xy2pix_nest(nside, ixp, iy, face_num)));
				res.set(6, new Long( xy2pix_nest(nside, ixp, iym, face_num)));
				res.set(7, new Long( xy2pix_nest(nside, ix, iym, face_num)));
				break;
			case 3: // NortWest edge
				other_face = 0 + ibm;
				ipo = (long) BitManipulation.MODULO(BitManipulation.invMSB( ipf), nsidesq); // NW-SE flip
				ixoiyo = pix2xy_nest(nside, ipo);
				ixo = ixoiyo[0];
				iyo = ixoiyo[1];
				res.set(2, new Long( xy2pix_nest(nside, ixo - 1, iyo,
						other_face)));
				res.set(3, new Long( (other_face * nsidesq + ipo)));
				res.set(4, new Long( xy2pix_nest(nside, ixo + 1, iyo,
						other_face)));
				res.set(0, new Long( xy2pix_nest(nside, ixm, iym, face_num)));
				res.set(1, new Long( xy2pix_nest(nside, ixm, iy, face_num)));
				res.set(5, new Long( xy2pix_nest(nside, ixp, iy, face_num)));
				res.set(6, new Long( xy2pix_nest(nside, ixp, iym, face_num)));
				res.set(7, new Long(xy2pix_nest(nside, ix, iym, face_num)));
				break;
			case 4: // SouthEast edge
				other_face = 8 + ib;
				ipo = (long) BitManipulation.MODULO(BitManipulation.invMSB( ipf), nsidesq); // SE-NW flip
				ixoiyo = pix2xy_nest(nside, ipo);
				ixo = ixoiyo[0];
				iyo = ixoiyo[1];
				res.set(0, new Long( xy2pix_nest(nside, ixo - 1, iyo,
						other_face)));
				res.set(1, new Long( xy2pix_nest(nside, ixm, iy, face_num)));
				res.set(2, new Long( xy2pix_nest(nside, ixm, iyp, face_num)));
				res.set(3, new Long( xy2pix_nest(nside, ix, iyp, face_num)));
				res.set(4, new Long( xy2pix_nest(nside, ixp, iyp, face_num)));
				res.set(5, new Long( xy2pix_nest(nside, ixp, iy, face_num)));
				res.set(6, new Long( xy2pix_nest(nside, ixo + 1, iyo,
						other_face)));
				res.set(7, new Long( (other_face * nsidesq + ipo)));
				break;
			case 5: // West corner
				other_face = 8 + ibm;
				arb_const = other_face * nsidesq + nsidesq - 1;
				res.set(0, new Long( (arb_const - 2)));
				res.set(1, new Long( arb_const));
				other_face = 4 + ibm;
				res.set(2, new Long( (other_face * nsidesq + local_magic1)));
				other_face = 0 + ibm;
				arb_const = other_face * nsidesq;
				res.set(3, new Long( arb_const));
				res.set(4, new Long( (arb_const + 1)));
				res.set(5, new Long( (ipix + 1)));
				res.set(6, new Long( (ipix - 1)));
				res.set(7, new Long( (ipix - 2)));
				break;
			case 6: //  North corner
				other_face = 0 + ibm;
				res.set(0, new Long( (ipix - 3)));
				res.set(1, new Long( (ipix - 1)));
				arb_const = other_face * nsidesq + local_magic1;
				res.set(2, new Long( (arb_const - 1)));
				res.set(3, new Long( arb_const));
				other_face = 0 + ib;
				arb_const = other_face * nsidesq + local_magic2;
				res.set(4, new Long( arb_const));
				res.set(5, new Long( (arb_const - 2)));
				res.set(6, new Long( (ipix - 2)));
				res.remove(7);
				break;
			case 7: // South corner
				other_face = 8 + ibm;
				arb_const = other_face * nsidesq + local_magic1;
				res.set(0, new Long( arb_const));
				res.set(1, new Long( (arb_const + 2)));
				res.set(2, new Long( (ipix + 2)));
				res.set(3, new Long( (ipix + 3)));
				res.set(4, new Long( (ipix + 1)));
				other_face = 8 + ib;
				arb_const = other_face * nsidesq + local_magic2;
				res.set(5, new Long( (arb_const + 1)));
				res.set(6, new Long( arb_const));
				res.remove(7);
				break;
			case 8: // East corner
				other_face = 8 + ib;
				arb_const = other_face * nsidesq + nsidesq - 1;
				res.set(0, new Long( (arb_const - 1)));
				res.set(1, new Long( (ipix - 1)));
				res.set(2, new Long( (ipix + 1)));
				res.set(3, new Long( (ipix + 2)));
				res.set(7, new Long( arb_const));
				other_face = 0 + ib;
				arb_const = other_face * nsidesq;
				res.set(4, new Long( (arb_const + 2)));
				res.set(5, new Long( arb_const));
				other_face = 4 + ibp;
				res.set(6, new Long( (other_face * nsidesq + local_magic2)));
				break;
			}
		} else { // South pole region
			switch (icase) {
			case 1: // North-East edge
				other_face = 4 + ibp;
				res.set(0, new Long( xy2pix_nest(nside, ixm, iym, face_num)));
				res.set(1, new Long( xy2pix_nest(nside, ixm, iy, face_num)));
				res.set(2, new Long( xy2pix_nest(nside, ixm, iyp, face_num)));
				res.set(3, new Long( xy2pix_nest(nside, ix, iyp, face_num)));
				res.set(7, new Long( xy2pix_nest(nside, ix, iym, face_num)));
				ipo = (long) BitManipulation.MODULO(BitManipulation.invLSB( ipf), nsidesq); // NE-SW flip
				ixoiyo = pix2xy_nest(nside, ipo);
				ixo = ixoiyo[0];
				iyo = ixoiyo[1];
				res.set(4, new Long( xy2pix_nest(nside, ixo, iyo + 1,
						other_face)));
				res.set(5, new Long( (other_face * nsidesq + ipo)));
				res.set(6, new Long( xy2pix_nest(nside, ixo, iyo - 1,
						other_face)));
				break;
			case 2: // SouthWest edge
				other_face = 8 + ibm;
				ipo = (long) BitManipulation.MODULO(BitManipulation.swapLSBMSB( ipf), nsidesq); // W-E flip
				ixoiyo = pix2xy_nest(nside, ipo);
				ixo = ixoiyo[0];
				iyo = ixoiyo[1];
				res.set(0, new Long( xy2pix_nest(nside, ixo - 1, iyo,
						other_face)));
				res.set(1, new Long( (other_face * nsidesq + ipo)));
				res.set(2, new Long( xy2pix_nest(nside, ixo + 1, iyo,
						other_face)));
				res.set(3, new Long( xy2pix_nest(nside, ix, iyp, face_num)));
				res.set(4, new Long( xy2pix_nest(nside, ixp, iyp, face_num)));
				res.set(5, new Long( xy2pix_nest(nside, ixp, iy, face_num)));
				res.set(6, new Long(xy2pix_nest(nside, ixp, iym, face_num)));
				res.set(7, new Long( xy2pix_nest(nside, ix, iym, face_num)));
				break;
			case 3: // NorthWest edge
				other_face = 4 + ib;
				ipo = (long) BitManipulation.MODULO(BitManipulation.invMSB( ipf), nsidesq);
				ixoiyo = pix2xy_nest(nside, ipo);
				ixo = ixoiyo[0];
				iyo = ixoiyo[1];
				res.set(0, new Long( xy2pix_nest(nside, ixm, iym, face_num)));
				res.set(1, new Long( xy2pix_nest(nside, ixm, iy, face_num)));
				res.set(2, new Long( xy2pix_nest(nside, ixo - 1, iyo,
						other_face)));
				res.set(3, new Long( (other_face * nsidesq + ipo)));
				res.set(4, new Long( xy2pix_nest(nside, ixo + 1, iyo,
						other_face)));
				res.set(5, new Long( xy2pix_nest(nside, ixp, iy, face_num)));
				res.set(6, new Long( xy2pix_nest(nside, ixp, iym, face_num)));
				res.set(7, new Long( xy2pix_nest(nside, ix, iym, face_num)));
				break;
			case 4: // SouthEast edge
				other_face = 8 + ibp;
				ipo = (long) BitManipulation.MODULO(BitManipulation.swapLSBMSB( ipf), nsidesq); // SE-NW
				// flip
				ixoiyo = pix2xy_nest(nside, ipo);
				ixo = ixoiyo[0];
				iyo = ixoiyo[1];
				res.set(0, new Long( xy2pix_nest(nside, ixo, iyo - 1,
						other_face)));
				res.set(1, new Long( xy2pix_nest(nside, ixm, iy, face_num)));
				res.set(2, new Long( xy2pix_nest(nside, ixm, iyp, face_num)));
				res.set(3, new Long( xy2pix_nest(nside, ix, iyp, face_num)));
				res.set(4, new Long( xy2pix_nest(nside, ixp, iyp, face_num)));
				res.set(5, new Long( xy2pix_nest(nside, ixp, iy, face_num)));
				res.set(6, new Long( xy2pix_nest(nside, ixo, iyo + 1,
						other_face)));
				res.set(7, new Long( (other_face * nsidesq + ipo)));
				break;
			case 5: // West corner
				other_face = 8 + ibm;
				arb_const = other_face * nsidesq + local_magic1;
				res.set(0, new Long( (arb_const - 2)));
				res.set(1, new Long( arb_const));
				other_face = 4 + ib;
				res.set(2, new Long( (other_face * nsidesq)));
				res.set(3, new Long( (other_face * nsidesq + 1)));
				res.set(4, new Long( (ipix + 1)));
				res.set(5, new Long( (ipix - 1)));
				res.set(6, new Long( (ipix - 2)));
				res.remove(7);
				break;
			case 6: //  North corner
				other_face = 4 + ib;
				res.set(0, new Long( (ipix - 3)));
				res.set(1, new Long((ipix - 1)));
				arb_const = other_face * nsidesq + local_magic1;
				res.set(2, new Long( (arb_const - 1)));
				res.set(3, new Long( arb_const));
				other_face = 0 + ib;
				res.set(4, new Long( (other_face * nsidesq)));
				other_face = 4 + ibp;
				arb_const = other_face * nsidesq + local_magic2;
				res.set(5, new Long( arb_const));
				res.set(6, new Long( (arb_const - 2)));
				res.set(7, new Long( (ipix - 2)));
				break;
			case 7: // South corner
				other_face = 8 + ib2;
				res.set(0, new Long( (other_face * nsidesq)));
				other_face = 8 + ibm;
				arb_const = other_face * nsidesq;
				res.set(1, new Long( arb_const));
				res.set(2, new Long( (arb_const + 1)));
				res.set(3, new Long( (ipix + 2)));
				res.set(4, new Long( (ipix + 3)));
				res.set(5, new Long( (ipix + 1)));
				other_face = 8 + ibp;
				arb_const = other_face * nsidesq;
				res.set(6, new Long( (arb_const + 2)));
				res.set(7, new Long( arb_const));
				break;
			case 8: // East corner
				other_face = 8 + ibp;
				res.set(1, new Long( (ipix - 1)));
				res.set(2, new Long( (ipix + 1)));
				res.set(3, new Long( (ipix + 2)));
				arb_const = other_face * nsidesq + local_magic2;
				res.set(6, new Long( arb_const));
				res.set(0, new Long( (arb_const - 2)));
				other_face = 4 + ibp;
				arb_const = other_face * nsidesq;
				res.set(4, new Long( (arb_const + 2)));
				res.set(5, new Long( arb_const));
				res.remove(7);
				break;
			}
		}
		return res;
	}

	/**
	 * returns the list of pixels in RING or NEST scheme with latitude in [phi0 -
	 * dpi, phi0 + dphi] on the ring iz in [1, 4*nside -1 ] The pixel id numbers
	 * are in [0, 12*nside^2 - 1] the indexing is in RING, unless nest is set to
	 * 1
	 * 
	 * @param nside 
	 *            long the map resolution
	 * @param iz 
	 *           long ring number
	 * @param phi0 
	 *            double
	 * @param dphi 
	 *            double
	 * @param nest 
	 *            boolean format flag
	 * @return ArrayList of pixels
	 * @throws IllegalArgumentException
	 * 
	 * Modified by N. Kuropatkin 07/09/2008  Corrected several bugs and make test of all cases.
	 * 
	 */
	public ArrayList InRing(long nside, long iz, double phi0, double dphi,
			boolean nest)  {
		boolean take_all = false;
		boolean to_top = false;
		boolean do_ring = true;
		boolean conservative = false;
//		String SID = "InRing:";
		double epsilon = Double.MIN_VALUE; // the constant to eliminate
		// java calculation jitter
		if (nest)
			do_ring = false;
		double shift = 0.;
		long ir = 0;
		long kshift, nr, ipix1, ipix2, nir1, nir2, ncap, npix;
		long ip_low = 0, ip_hi = 0, in, inext, nir;
		ArrayList res = new ArrayList();
		npix = 12 * nside * nside; // total number of pixels
		ncap = 2 * nside * (nside - 1); // number of pixels in the north polar
										// cap
		double phi_low = BitManipulation.MODULO(phi0 - dphi, TWOPI) - epsilon; // phi min,
																  // excluding
																  // 2pi period
		double phi_hi = BitManipulation.MODULO(phi0 + dphi, TWOPI) + epsilon;

//
		if (Math.abs(dphi - PI) < epsilon)  take_all = true;

		/* identifies ring number */
		if ((iz >= nside) && (iz <= 3 * nside)) { // equatorial region
			ir = iz - nside + 1; // in [1, 2*nside + 1]
			ipix1 = ncap + 4 * nside * (ir - 1); // lowest pixel number in the
											     // ring
			ipix2 = ipix1 + 4 * nside - 1; // highest pixel number in the ring
			kshift = (long) BitManipulation.MODULO(ir, 2.);

			nr = nside * 4;
		} else {
			if (iz < nside) { // north pole
				ir = iz;
				ipix1 = 2 * ir * (ir - 1); // lowest pixel number
				ipix2 = ipix1 + 4 * ir - 1; // highest pixel number
			} else { // south pole
				ir = 4 * nside - iz;

				ipix1 = npix - 2 * ir * (ir + 1); // lowest pixel number
				ipix2 = ipix1 + 4 * ir - 1;       // highest pixel number
			}
			nr = ir * 4;
			kshift = 1;
		}
		
		// Construct the pixel list
		if (take_all) {
			nir = ipix2 - ipix1 + 1;
			if (do_ring) {
				long ind = 0;
				for (long i =  ipix1; i <= ipix2; i++) {
					res.add((int) ind, new Long(i));
					ind++;
				}
			} else {
//				in = ring2nest(nside, ipix1);
//				res.add(0, new Long( in));
				int ind = 0;
				for (int i = 0; i < nir; i++) {
//					inext = next_in_line_nest(nside, in);
//					in = inext;
					in = ring2nest(nside, ipix1 + i);
					res.add( i, new Long(in));

				}
			}

			return res;
		}
		
		shift = kshift / 2.0;

		// conservative : include every intersected pixel, even if the
		// pixel center is out of the [phi_low, phi_hi] region
		if (conservative) {
			ip_low = (long) Math.round((nr * phi_low) / TWOPI - shift);
			ip_hi = (long) Math.round((nr * phi_hi) / TWOPI - shift);

			ip_low = (long) BitManipulation.MODULO(ip_low, nr); // in [0, nr - 1]
			ip_hi = (long) BitManipulation.MODULO(ip_hi, nr); // in [0, nr - 1]
//			System.out.println("ip_low="+ip_low+" ip_hi="+ip_hi);
		} else { // strict: includes only pixels whose center is in
			//                                                    [phi_low,phi_hi]

			ip_low = (long) Math.ceil((nr * phi_low) / TWOPI - shift);
			ip_hi = (long) Math.floor((nr * phi_hi) / TWOPI - shift);
			if (ip_low == ip_hi + 1)
				ip_low = ip_hi;

			if ((ip_low - ip_hi == 1) && (dphi * nr < PI)) {
				// the interval is too small ( and away from pixel center)
				// so no pixels is included in the list
				System.out
						.println("the interval is too small and avay from center");
				return res; // return empty list
			}

			ip_low = Math.min(ip_low, nr - 1);
			ip_hi = Math.max(ip_hi, 0);
		}

		//
		if (ip_low > ip_hi)
			to_top = true;

		if (to_top) {
			ip_low += ipix1;
			ip_hi += ipix1;
			nir1 = ipix2 - ip_low + 1;

			nir2 = ip_hi + 1;
			nir = nir1 + nir2;
			if (do_ring) {
				int ind = 0;
				for (long i =  ip_low; i <= ipix2; i++) {
					res.add(ind, new Long(i));
					ind++;
				}
				//				ind = nir1;
				for (long i =  ipix1; i <= ip_hi; i++) {
					res.add(ind, new Long(i));
					ind++;
				}
			} else {
				int ind = 0;
				for (long i = ip_low; i <= ipix2; i++) {
					in = ring2nest(nside, i);
					res.add(ind, new Long(in));
					ind++;
				}
				for (long i =  ipix1; i <= ip_hi; i++) {
					in = ring2nest(nside, i);
					res.add(ind, new Long(in));
					ind++;
				}
				 
			}
		} else {
			nir = ip_hi - ip_low + 1;
//			System.out.println("nir="+nir+" ip_low="+ip_low+" ip_hi="+ip_hi+" ipix1="+
//					ipix1+" ipix2="+ipix2);
			//
			// Special case when ip_low < 0
			//
			if (ip_low < 0 ){
				ip_low = Math.abs(ip_low) ;
				nir1 = ip_low;
				nir2 = ip_hi + 1;
				nir = nir1 + nir2;
				if (do_ring) {
					int ind = 0;
					for (long i =  0; i < ip_low; i++) {
						res.add(ind, new Long(ipix2 - i));
						ind++;
					}
					for (long i =  0; i <= ip_hi; i++) {
						res.add(ind, new Long(ipix1 + i));
						ind++;
					}
					
				} else {
					
					int ind = 0;
					for (int i = 0; i < ip_low; i++) {
						in = ring2nest(nside, ipix2 - i);						
						res.add(ind, new Long(in));
						ind++;
					}
					for (long i =  0; i <= ip_hi; i++) {
						in = ring2nest(nside, ipix1 + i);
						res.add(ind, new Long(in));
						ind++;
					}

				}
//				System.out.println("nir="+nir+" ip_low="+ip_low+" ip_hi="+ip_hi+" ipix1="+
//						ipix1+" ipix2="+ipix2);
//				for (int i=0; i< res.size(); i++) {
//					System.out.print(" "+ res.get(i));
//				}
//				System.out.println();
				return res;
				
			}
			ip_low += ipix1;
			ip_hi += ipix1;
			if (do_ring) {
				int ind = 0;
				for (long i =  ip_low; i <= ip_hi; i++) {
					res.add(ind, new Long(i));
					ind++;
				}
			} else {
//				in = ring2nest(nside, ip_low);
//				res.add(0, new Long(in));
				int ind = 0;
				for (long i = ip_low; i <= ip_hi; i++) {
//					inext = next_in_line_nest(nside, in);
//					in = inext;
					in = ring2nest(nside, i);
					res.add(ind, new Long(in));
					ind++;
				}
			}

		}

		return res;
	}

	/**
	 * calculates the pixel that lies on the East side (and the same
	 * latitude) as the given NESTED pixel number - ipix
	 * 
	 * @param nside 
	 *            long resolution
	 * @param ipix 
	 *            long pixel number
	 * @return  long next pixel in line
	 * @throws IllegalArgumentException
	 */
	public long next_in_line_nest(long nside, long ipix)  {
		long npix, ipf, ipo, ix, ixp, iy, iym, ixo, iyo, face_num, other_face;
		long ia, ib, ibp, nsidesq;
		long ibm, ib2;
        int icase;
		long local_magic1, local_magic2;
		long[] ixiy = new long[2];
		long inext = 0; // next in line pixel in Nest scheme
		String SID = "next_in_line:";
		if ((nside < 1) || (nside > ns_max)) {
			throw new IllegalArgumentException(SID + " nside should be power of 2 >0 and < "+ns_max);
		}
		nsidesq = nside * nside;
		npix = 12 * nsidesq; // total number of pixels
		if ((ipix < 0) || (ipix > npix - 1)) {
			throw new IllegalArgumentException(SID + " ipix out of range defined by nside");
		}
		// initiates array for (x,y) -> pixel number -> (x,y) mapping
		if (x2pix[xmax-1] <= 0)
			mk_xy2pix();
		local_magic1 = (nsidesq - 1) / 3;
		local_magic2 = 2 * local_magic1;
		face_num = ipix / nsidesq;
		ipf = (long) BitManipulation.MODULO(ipix, nsidesq); // Pixel number in face
		ixiy = pix2xy_nest(nside, ipf);
		ix = ixiy[0];
		iy = ixiy[1];
		ixp = ix + 1;
		iym = iy - 1;
		boolean sel = false;
		icase = -1; // iside the nest flag
		// Exclude corners
		if (ipf == local_magic2) { // west coirner
			inext = ipix - 1;
			return inext;
		}
		if ((ipf == nsidesq - 1) && !sel) { // North corner
			icase = 6;
			sel = true;
		}
		if ((ipf == 0) && !sel) { // Siuth corner
			icase = 7;
			sel = true;
		}
		if ((ipf == local_magic1) && !sel) { // East corner
			icase = 8;
			sel = true;
		}
		// Detect edges
		if (((ipf & local_magic1) == local_magic1) && !sel) { // North-East
			icase = 1;
			sel = true;
		}
		if (((ipf & local_magic2) == 0) && !sel) { // South-East
			icase = 4;
			sel = true;
		}
		if (!sel) { // iside a face
			inext = xy2pix_nest(nside, ixp, iym, face_num);
			return inext;
		}
		//
		ia = face_num / 4; // in [0,2]
		ib = (long) BitManipulation.MODULO(face_num, 4); // in [0,3]
		ibp = (long) BitManipulation.MODULO(ib + 1, 4);
		ibm = (long) BitManipulation.MODULO(ib + 4 - 1, 4);
		ib2 = (long) BitManipulation.MODULO(ib + 2, 4);

		if (ia == 0) { // North pole region
			switch (icase) {
			case 1:
				other_face = 0 + ibp;
				ipo = (long) BitManipulation.MODULO(BitManipulation.swapLSBMSB( ipf), nsidesq);
				inext = other_face * nsidesq + ipo;
				break;
			case 4:
				other_face = 4 + ibp;
				ipo = (long) BitManipulation.MODULO(BitManipulation.invMSB( ipf), nsidesq); // SE-NW flip

				ixiy = pix2xy_nest(nside, ipo);
				ixo = ixiy[0];
				iyo = ixiy[1];

				inext = xy2pix_nest(nside, ixo + 1, iyo, other_face);

				break;
			case 6: // North corner
				other_face = 0 + ibp;
				inext = other_face * nsidesq + nsidesq - 1;
				break;
			case 7:
				other_face = 4 + ibp;
				inext = other_face * nsidesq + local_magic2 + 1;
				break;
			case 8:
				other_face = 0 + ibp;
				inext = other_face * nsidesq + local_magic2;
				break;
			}

		} else if (ia == 1) { // Equatorial region
			switch (icase) {
			case 1: // NorthEast edge
				other_face = 0 + ib;
//                System.out.println("ipf="+ipf+" nsidesq="+nsidesq+" invLSB="+BitManipulation.invLSB(ipf));
				ipo = (long) BitManipulation.MODULO((double)BitManipulation.invLSB( ipf), (double)nsidesq); // NE-SW flip
//                System.out.println(" ipo="+ipo);
                
				ixiy = pix2xy_nest(nside, ipo);
				ixo = ixiy[0];
				iyo = ixiy[1];
				inext = xy2pix_nest(nside, ixo, iyo - 1, other_face);
				break;
			case 4: // SouthEast edge
				other_face = 8 + ib;
				ipo = (long) BitManipulation.MODULO(BitManipulation.invMSB( ipf), nsidesq);
				ixiy = pix2xy_nest(nside, ipo);
				inext = xy2pix_nest(nside, ixiy[0] + 1, ixiy[1], other_face);
				break;
			case 6: // Northy corner
				other_face = 0 + ib;
				inext = other_face * nsidesq + local_magic2 - 2;
				break;
			case 7: // South corner
				other_face = 8 + ib;
				inext = other_face * nsidesq + local_magic2 + 1;
				break;
			case 8: // East corner
				other_face = 4 + ibp;
				inext = other_face * nsidesq + local_magic2;
				break;

			}
		} else { // South pole region
			switch (icase) {
			case 1: // NorthEast edge
				other_face = 4 + ibp;
				ipo = (long) BitManipulation.MODULO(BitManipulation.invLSB( ipf), nsidesq); // NE-SW flip
				ixiy = pix2xy_nest(nside, ipo);
				inext = xy2pix_nest(nside, ixiy[0], ixiy[1] - 1, other_face);
				break;
			case 4: // SouthEast edge
				other_face = 8 + ibp;
				ipo = (long) BitManipulation.MODULO(BitManipulation.swapLSBMSB( ipf), nsidesq); // E-W flip
				inext = other_face * nsidesq + ipo; // (8)
				break;
			case 6: // North corner
				other_face = 4 + ibp;
				inext = other_face * nsidesq + local_magic2 - 2;
				break;
			case 7: // South corner
				other_face = 8 + ibp;
				inext = other_face * nsidesq;
				break;
			case 8: // East corner
				other_face = 8 + ibp;
				inext = other_face * nsidesq + local_magic2;
				break;
			}
		}
		return inext;
	}

	/**
	 * renders the pixel number pix (NESTED scheme) for a pixel which contains a
	 * point on a sphere at coordinate vector (x,y,z), given the map resolution
	 * parameter nside.
	 * 
	 * The computation is made to the highest resolution available (nside=ns_max)
	 * and then degraded to that requared (by Integer division) this doesn't
	 * cost much, and it makes sure that the treatment of round-off will be
	 * consistent for every resolution
	 * 
	 * @param nside
	 *            long the map resolution
	 * @param vector
	 *            Vewctor3d the input vector
	 * @return pixel long
	 * @throws IllegalArgumentException
	 */
	public long vect2pix_nest(long nside, Vector3d vector)  {
		long pix = 0;
		double z, za, tt, tp, tmp, dnorm, phi;
		long jp, jm, ifp, ifm, face_num, ix, iy, ix_low, ix_hi;
		long iy_low, iy_hi, ipf, ntt;
		String SID = " vect2pix_nest:";
		/*                      */
		if (nside < 1 || nside > ns_max) {
			throw new IllegalArgumentException(SID + " nside should be power of 2 >0 and < ns_max");
		}
		if (x2pix[xmax-1] <= 0)
			mk_xy2pix();

		dnorm = vector.length();
		z = vector.z / dnorm;
		phi = 0.; // phi in [-pi,pi]
		if (vector.x != 0. || vector.y != 0.)
			phi = Math.atan2(vector.y, vector.x);
		za = Math.abs(z);
		if (phi < 0.)
			phi += TWOPI; // phi in [0,2pi]
		tt = 2. * phi / PI; // in [0,4]
		if (za <= twothird) { // Equatorial region
			/*
			 * the index of edge lines increases when the longitude = phi goes
			 * up
			 */
			jp = (long) (ns_max * (0.5 + tt - z * 0.75)); // ascending edge line
			// index
			jm = (long) (ns_max * (0.5 + tt + z * 0.75)); // descending edge line
			// index
			/* find the face */
			ifp = jp / ns_max; // in [0,4]
			ifm = jm / ns_max;
			if (ifp == ifm) { // faces 4 to 7
				face_num = (long) BitManipulation.MODULO(ifp, 4) + 4;
			} else if (ifp < ifm) { // (half-) faces 0 to 3
				face_num = (long) BitManipulation.MODULO(ifp, 4);
			} else { // (half-) faces 8 to 11
				face_num = (long) BitManipulation.MODULO(ifm, 4) + 8;
			}

			ix = (long) BitManipulation.MODULO(jm, ns_max);
			iy = ns_max - (long) BitManipulation.MODULO(jp, ns_max) - 1;
		} else { // polar region, za > 2/3
			ntt = (long) tt;
			if (ntt >= 4)
				ntt = 3;
			tp = tt - ntt;
			tmp = Math.sqrt(3.0 * (1.0 - za)); // in [0,1]
			/*
			 * the index of edge lines increases when distance from the closest
			 * pole goes up
			 */
			jp = (long) (ns_max * tp * tmp); // line going toward the pole has
			// phi increases
			jm = (long) (ns_max * (1.0 - tp) * tmp); // that one goes away of the
			// closest pole
			jp = Math.min(ns_max - 1, jp); // for points too close to the
			// boundary
			jm = Math.min(ns_max - 1, jm);
			/* finds the face and pixel's (x,y) */
			if (z >= 0) {
				face_num = ntt; // in [0,3]
				ix = ns_max - jm - 1;
				iy = ns_max - jp - 1;
			} else {
				face_num = ntt + 8; // in [8,11]
				ix = jp;
				iy = jm;
			}
		}
		ix_low = (long) BitManipulation.MODULO(ix, xmax);
		ix_hi = ix / xmax;
		iy_low = (long) BitManipulation.MODULO(iy, xmax);
		iy_hi = iy / xmax;

		ipf = (x2pix[(int) (ix_hi + 1)] + y2pix[(int) (iy_hi + 1)]) * (xmax * xmax)
				+ (x2pix[(int) (ix_low + 1)] + y2pix[(int) (iy_low + 1)]);
		ipf = ipf / ((ns_max / nside) * (ns_max / nside)); // in [0, nside**2
		// -1]
		pix = ipf + face_num * nside * nside; // in [0, 12*nside**2 -1]
		return pix;
	}

	/**
	 * gives the pixel number ipix (NESTED) corresponding to ix, iy and face_num
	 * 
	 * @param nside 
	 *            the map resolution parameter
	 * @param ix 
	 *            Integer x coordinate
	 * @param iy 
	 *            Integer y coordinate
	 * @param face_num 
	 *            long face number
	 * @return  long pixel number ipix
	 * @throws IllegalArgumentException
	 */
	private long xy2pix_nest(long nside, long ix, long iy, long face_num){
		long res = 0;
		long ix_low, ix_hi, iy_low, iy_hi, ipf;
		String SID = "xy2pix_nest:";
		//
		if ((nside < 1) || (nside > ns_max)) {
			throw new IllegalArgumentException(SID + " nside should be power of 2 >0 and < "+ns_max);
		}
//		if ((ix < 0) || (ix > nside - 1)) {
//			throw new IllegalArgumentException(SID + " ix out of range [0, nside-1]");
//		}
//		if ((iy < 0) || (iy > nside - 1)) {
//			throw new IllegalArgumentException(SID + " iy out of range [0, nside-1]");
//		}
		if (x2pix[xmax-1] <= 0)
			mk_xy2pix();
		ix_low = (long) BitManipulation.MODULO(ix, xmax);
		ix_hi = ix / xmax;
		iy_low = (long) BitManipulation.MODULO(iy, xmax);
		iy_hi = iy / xmax;

		ipf = (x2pix[(int) (ix_hi + 1)] + y2pix[(int) (iy_hi + 1)]) * xmax * xmax
				+ (x2pix[(int) (ix_low + 1)] + y2pix[(int) (iy_low + 1)]);
		res = ipf + face_num * nside * nside; // in [0, 12*nside^2 - 1]

		return res;
	}

	/**
	 * gives the x,y coordinates in a face from pixel number within the face
	 * (NESTED) schema.
	 * 
	 * @param nside 
	 *            long resolution parameter
	 * @param ipf 
	 *            long pixel number
	 * @return ixiy  long[] contains x and y coordinates
	 * @throws IllegalArgumentException
	 */
	private long[] pix2xy_nest(long nside, long ipf)  {
		long[] ixiy = { 0, 0 };
		long ip_low, ip_trunc, ip_med, ip_hi;
		String SID = "pix2xy_nest:";
//        System.out.println(" ipf="+ipf+" nside="+nside);
		if ((nside < 1) || (nside > ns_max)) {
			throw new IllegalArgumentException(SID + " nside should be power of 2 >0 and < "+ns_max);
		}
		if ((ipf < 0) || (ipf > nside * nside - 1)) {
			throw new IllegalArgumentException(SID + " ipix out of range defined by nside");
		}
		if (pix2x[pixmax] <= 0)
			mk_pix2xy();
		ip_low = (long) BitManipulation.MODULO(ipf, pixmax); // contents of last 15 bits
		ip_trunc = ipf / pixmax; // truncation of the last 15 bits
		ip_med = (long) BitManipulation.MODULO(ip_trunc, pixmax); // next 15 bits
		ip_hi = ip_trunc / pixmax; // select high 15 bits

		long ix = pixmax * pix2x[(int) ip_hi] + xmid * pix2x[(int) ip_med] + pix2x[(int) ip_low];
		long iy = pixmax * pix2y[(int) ip_hi] + xmid * pix2y[(int) ip_med] + pix2y[(int) ip_low];
		ixiy[0] = ix;
		ixiy[1] = iy;
		return ixiy;
	}

	/**
	 * creates an array of pixel numbers pix2x from x and y coordinates in the
	 * face. Suppose NESTED scheme of pixel ordering Bits corresponding to x and
	 * y are interleaved in the pixel number in even and odd bits.
	 */
	private void mk_pix2xy() {
		long kpix, jpix, ix, iy, ip, id;
		boolean flag = true;
		for (kpix = 0; kpix <= pixmax; kpix++) { // loop on pixel numbers
			jpix = kpix;
			ix = 0;
			iy = 0;
			ip = 1; // bit position in x and y

			while (jpix != 0) { // go through all the bits

				id = (long) BitManipulation.MODULO(jpix, 2); // bit value (in kpix), goes in x
				jpix /= 2;
				ix += id * ip;

				id = (long) BitManipulation.MODULO(jpix, 2); // bit value, goes in iy
				jpix /= 2;
				iy += id * ip;

				ip *= 2; // next bit (in x and y )
			}
 
			pix2x[(int) kpix] = ix; // in [0,pixmax]
			pix2y[(int) kpix] = iy; // in [0,pixmax]
			

		}
	}

	/**
	 * converts pixel number from ring numbering schema to the nested one
	 * 
	 * @param nside 
	 *            long resolution
	 * @param ipring long pixel number in ring schema
	 * @return long pixel number in nest schema
	 * @throws IllegalArgumentException
	 */
	public long ring2nest(long nside, long ipring)  {
		long ipnest = 0;
		double fihip;
		double hip;
		long npix, nl2, nl4, ncap, ip, iphi, ipt, ipring1, kshift, face_num;
		long nr, irn, ire, irm, irs, irt, ifm, ifp, ix, iy, ix_low, ix_hi, iy_low;
		long iy_hi, ipf;
		long[] jrll = { 0, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 }; // in units of
		// nside
		long[] jpll = { 0, 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 }; // in units of
		// nside/2
		String SID = "ring2nest:";
		//
		face_num = 0;
		if ((nside < 1) || (nside > ns_max)) {
			throw new IllegalArgumentException(SID + " nside should be power of 2 >0 and < "+ns_max);
		}
		npix = 12 * nside * nside; // total number of points

		if ((ipring < 0) || (ipring > npix - 1)) {
			throw new IllegalArgumentException(SID + " ipring out of range [0,npix-1]");
		}
		if (x2pix[xmax-1] <= 0)
			mk_xy2pix();

		nl2 = 2 * nside;
		nl4 = 4 * nside;
		ncap = nl2 * (nside - 1); // points in each polar cap, =0 for nside = 1
		ipring1 = ipring + 1;
		// finds the ring number, the position of the ring and the face number
		if (ipring1 <= ncap) { // north polar cap
			hip = ipring1 / 2.0;
			fihip = Math.floor(hip);
			irn = (long)Math.floor( Math.sqrt(hip - Math.sqrt(fihip))) + 1; // counted from
			// north pole
			iphi = ipring1 - 2 * irn * (irn - 1);

			kshift = 0;
			nr = irn; // 1/4 of the number of points on the current ring
			face_num = (iphi - 1) / irn; // in [0,3 ]
			
		} else if (ipring1 <= nl2 * (5 * nside + 1)) { // equatorial region		
			ip = ipring1 - ncap - 1;
			irn = (long)Math.floor(ip / nl4) + nside; // counted from north pole
			iphi = (long) BitManipulation.MODULO(ip, nl4) + 1;

			kshift = (long) BitManipulation.MODULO(irn + nside, 2); // 1 if odd 0
			// otherwise
			nr = nside;
			ire = irn - nside + 1; // in [1, 2*nside+1]
			irm = nl2 + 2 - ire;
			ifm = (iphi - ire / 2 + nside - 1) / nside; // face boundary
			ifp = (iphi - irm / 2 + nside - 1) / nside;
			if (ifp == ifm) {
				face_num = (long) BitManipulation.MODULO(ifp, 4.) + 4;
			} else if (ifp + 1 == ifm) { // (half-) faces 0 to 3
				face_num = ifp;
			} else if (ifp - 1 == ifm) { // (half-) faces 8 to 11
				face_num = ifp + 7;
			}
		
		
		} else { // south polar cap
			
			ip = npix - ipring1 + 1;
			hip = ip / 2.0;
			fihip = Math.floor(hip);
			irs = (long)Math.floor( Math.sqrt(hip - Math.sqrt(fihip))) + 1;
			iphi = 4 * irs + 1 - (ip - 2 * irs * (irs - 1));
			kshift = 0;
			nr = irs;
			irn = nl4 - irs;
			face_num = (iphi - 1) / irs + 8; // in [8,11]
			
			
		}
		// finds the (x,y) on the face
		
		
//
		irt = irn - jrll[(int) (face_num + 1)] * nside + 1; // in [-nside+1,0]
		ipt = 2 * iphi - jpll[(int) (face_num + 1)] * nr - kshift - 1; // in [-nside+1,
		// nside-1]
//
		if (ipt >= nl2){
			ipt = ipt - 8*nside; // for the face #4		
		}
		ix = (ipt - irt) / 2;
		iy = -(ipt + irt) / 2;

		ix_low = (long) BitManipulation.MODULO(ix, xmax);
		ix_hi = ix / xmax;
		iy_low = (long) BitManipulation.MODULO(iy, xmax);
		iy_hi = iy / xmax;

          //
		
		ipf = (x2pix[(int) (ix_hi + 1)] + y2pix[(int) (iy_hi + 1)]) * xmax * xmax
				+ (x2pix[(int) (ix_low + 1)] + y2pix[(int) (iy_low + 1)]); // in [0, nside**2 -1]
		ipnest = ipf + face_num * nside * nside; // in [0, 12*nside**2 -1]
		
		return ipnest;

	}

	/**
	 * converts from NESTED to RING pixel numbering
	 * 
	 * @param nside 
	 *            long resolution
	 * @param ipnest
	 *            long NEST pixel number
	 * @return ipring  long RING pixel number
	 * @throws IllegalArgumentException
	 */
	public long nest2ring(long nside, long ipnest)  {
		long res = 0;
		long npix, npface, face_num, ncap, n_before, ipf, ip_low, ip_trunc;
		long ip_med, ip_hi, ix, iy, jrt, jr, nr, jpt, jp, kshift, nl4;
//		long[] ixiy = { 0, 0 };
		// coordinates of lowest corner of each face
		long[] jrll = { 0, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 }; // in units of
		// nside
		long[] jpll = { 0, 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 }; // in units of
		// nside/2
		String SID = "nest2ring:";
		//
		if ((nside < 1) || (nside > ns_max)) {
			throw new IllegalArgumentException(SID + " nside should be power of 2 >0 and < ns_max");
		}
		npix = 12 * nside * nside;
		if ((ipnest < 0) || (ipnest > npix - 1)) {
			throw new IllegalArgumentException(SID + " ipnest out of range [0,npix-1]");
		}
		if (pix2x[pixmax-1] <= 0)
			mk_pix2xy();
		ncap = 2 * nside * (nside - 1); // number of points in the North polar
		// cap
		nl4 = 4 * nside;
		// finds the face and the number in the face
		npface = nside * nside;

		face_num = ipnest / npface; // face number in [0,11]
		if (ipnest >= npface) {
			ipf = (long) BitManipulation.MODULO(ipnest, npface); // pixel number in the face
		} else {
			ipf = ipnest;
		}

		// finds the x,y on the face
		//  from the pixel number
		ip_low = (long) BitManipulation.MODULO(ipf, pixmax); // last 15 bits
		if (ip_low < 0)
			ip_low = -ip_low;

		ip_trunc = ipf / pixmax; // truncate last 15 bits
		ip_med = (long) BitManipulation.MODULO(ip_trunc, pixmax); // next 15 bits
		if (ip_med < 0)
			ip_med = -ip_med;
		ip_hi = ip_trunc / pixmax; // high 15 bits

		ix = pixmax * pix2x[(int) ip_hi] + xmid * pix2x[(int) ip_med] + pix2x[(int) ip_low];
		iy = pixmax * pix2y[(int) ip_hi] + xmid * pix2y[(int) ip_med] + pix2y[(int) ip_low];

		// transform this in (horizontal, vertical) coordinates
		jrt = ix + iy; // vertical in [0,2*(nside -1)]
		jpt = ix - iy; // horizontal in [-nside+1, nside - 1]
		// calculate the z coordinate on the sphere
		jr = jrll[(int) (face_num + 1)] * nside - jrt - 1; // ring number in [1,4*nside
		// -1]
		nr = nside; // equatorial region (the most frequent)
		n_before = ncap + nl4 * (jr - nside);
		kshift = (long) BitManipulation.MODULO(jr - nside, 2);
		if (jr < nside) { // north pole region
			nr = jr;
			n_before = 2 * nr * (nr - 1);
			kshift = 0;
		} else if (jr > 3 * nside) { // south pole region
			nr = nl4 - jr;
			n_before = npix - 2 * (nr + 1) * nr;
			kshift = 0;
		}
		// computes the phi coordinate on the sphere in [0,2pi]
		jp = (jpll[(int) (face_num + 1)] * nr + jpt + 1 + kshift) / 2; // 'phi' number
		// in ring
		// [1,4*nr]
		if (jp > nl4)
			jp -= nl4;
		if (jp < 1)
			jp += nl4;
		res = n_before + jp - 1; // in [0, npix-1]
		return res;
	}

	/**
	 * fills arrays x2pix and y2pix giving the number of the pixel laying in
	 * (x,y). x and y are in [1,512] the pixel number is in [0, 512**2 -1]
	 * 
	 * if i-1 = sum_p=0 b_p*2^p then ix = sum+p=0 b_p*4^p iy = 2*ix ix + iy in
	 * [0,512**2 -1]
	 * 
	 */
	private static void mk_xy2pix() {
		long k, ip, id;

		for (int i = 1; i <= xmax; i++) {
			long j = i - 1;
			k = 0;
			ip = 1;
			while (j != 0) {
				id = (long) BitManipulation.MODULO(j, 2);
				j /= 2;
				k += ip * id;
				ip *= 4;
			}
			x2pix[i] = k;
			y2pix[i] = 2 * k;
			
		}

	}

	/**
	 * returns the ring number in {1, 4*nside - 1} calculated from z coordinate
	 * 
	 * @param nside 
	 *            long resolution
	 * @param z 
	 *            double z coordinate
	 * @return long ring number
	 */
	public long RingNum(long nside, double z) {
		long iring = 0;
		/* equatorial region */

		iring = (long) Math.round(nside * (2.0 - 1.5 * z));
		/* north cap */
		if (z > twothird) {
			iring = (long) Math.round(nside * Math.sqrt(3.0 * (1.0 - z)));
			if (iring == 0)
				iring = 1;
		}
		/* south cap */
		if (z < -twothird) {
			iring = (long) Math.round(nside * Math.sqrt(3.0 * (1.0 + z)));
			if (iring == 0)
				iring = 1;
			iring = 4 * nside - iring;
		}
		return iring;
	}

	/**
	 * calculates vector corresponding to angles theta (co-latitude
	 * measured from North pole, in [0,pi] radians) phi (longitude measured
	 * eastward in [0,2pi] radians) North pole is (x,y,z) = (0, 0, 1)
	 * 
	 * @param theta double
	 * @param phi double
	 * @return Vector3d
	 * @throws IllegalArgumentException
	 */
	public Vector3d Ang2Vec(double theta, double phi)  {
		double PI = Math.PI;
		String SID = "Ang2Vec:";
		Vector3d v;
		if ((theta < 0.0) || (theta > PI)) {
			throw new IllegalArgumentException(SID + " theta out of range [0.,PI]");
		}
		double stheta = Math.sin(theta);
		double x = stheta * Math.cos(phi);
		double y = stheta * Math.sin(phi);
		double z = Math.cos(theta);
		v = new Vector3d(x, y, z);
		return v;
	}

	/**
	 * converts a Vector3d in a tuple of angles tup[0] = theta 
	 * co-latitude measured from North pole, in [0,PI] radians, tup[1] = phi 
	 * longitude measured eastward, in [0,2PI] radians
	 * 
	 * @param v
	 *            Vector3d
	 * @return double[] out_tup out_tup[0] = theta out_tup[1] = phi
	 */
	public double[] Vect2Ang(Vector3d v) {
		double[] out_tup = new double[2];
		double norm = v.length();
		double z = v.z / norm;
		double theta = Math.acos(z);
		double phi = 0.;
		if ((v.x != 0.) || (v.y != 0)) {
			phi = Math.atan2(v.y, v.x); // phi in [-pi,pi]
		}
		if (phi < 0)
			phi += 2.0 * Math.PI; // phi in [0, 2pi]
//		phi += Math.PI;
		out_tup[0] = theta;
		out_tup[1] = phi;
		return out_tup;
	}

	/**
	 * returns nside such that npix = 12*nside^2,  nside should be
	 * power of 2 and smaller than ns_max if not return -1
	 * 
	 * @param npix
	 *            long the number of pixels in the map
	 * @return long nside the map resolution parameter
	 */
	public long Npix2Nside(long npix) {
		long nside = 0;
		long npixmax = 12 *(long) ns_max *(long) ns_max;
 
		String SID = "Npix2Nside:";
		nside = (long) Math.rint(Math.sqrt(npix / 12));
		if (npix < 12) {
			throw new IllegalArgumentException(SID + " npix is too small should be > 12");
		}
		if (npix > npixmax) {
			throw new IllegalArgumentException(SID + " npix is too large > 12 * ns_max^2");
		}
		double fnpix = 12.0 * nside * nside;
		if (Math.abs(fnpix - npix) > 1.0e-2) {
			throw new IllegalArgumentException(SID + "  npix is not 12*nside*nside");
		}
		double flog = Math.log((double) nside) / Math.log(2.0);
		double ilog = Math.rint(flog);
		if (Math.abs(flog - ilog) > 1.0e-6) {
			throw new IllegalArgumentException(SID + "  nside is not power of 2");
		}
		return nside;
	}

	/**
	 * calculates npix such that npix = 12*nside^2 ,nside should be
	 * a power of 2, and smaller than ns_max otherwise return -1 
	 * 
	 * @param nside
	 *            long the map resolution
	 * @return npix long the number of pixels in the map
	 */
	public long Nside2Npix(long nside) {

		long[] nsidelist = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048,
				4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576,
                2097152, 4194304};

		long res = 0;
		String SID = "Nside2Npix:";
		if (Arrays.binarySearch(nsidelist, nside) < 0) {
			throw new IllegalArgumentException(SID + " nside should be >0, power of 2, <"+ns_max);
		}
		res = 12 * nside * nside;
		return res;
	}

	/**
	 * calculates the surface of spherical triangle defined by
	 * vertices v1,v2,v3 Algorithm: finds triangle sides and uses l'Huilier
	 * formula to compute "spherical excess" = surface area of triangle on a
	 * sphere of radius one see, eg Bronshtein, Semendyayev Eq 2.86 half
	 * perimeter hp = 0.5*(side1+side2+side3) l'Huilier formula x0 = tan( hp/2.)
	 * x1 = tan((hp - side1)/2.) x2 = tan((hp - side2)/2.) x3 = tan((hp -
	 * side3)/2.)
	 * 
	 * @param v1 Vector3d
	 * @param v2 Vector3d
	 * @param v3 Vector3d vertices of the triangle
	 * @return  double the triangle surface
	 * @throws Exception
	 *  
	 */
	public double SurfaceTriangle(Vector3d v1, Vector3d v2, Vector3d v3)
			throws Exception {
		double res = 0.;
		double side1 = AngDist(v2, v3) / 4.0;
		double side2 = AngDist(v3, v1) / 4.0;
		double side3 = AngDist(v1, v2) / 4.0;
		double x0 = Math.tan(side1 + side2 + side3);
		double x1 = Math.tan(side2 + side3 - side1);
		double x2 = Math.tan(side1 + side3 - side2);
		double x3 = Math.tan(side1 + side2 - side3);
		res = 4.0 * Math.atan(Math.sqrt(x0 * x1 * x2 * x3));

		return res;
	}

	/**
	 * calculates angular distance (in radians) between 2 Vectors
	 * v1 and v2 In general dist = acos(v1.v2) except if the vectors are almost
	 * aligned
	 * 
	 * @param v1 Vector3d
	 * @param v2 Vector3d
	 * @return double dist in radians
	 * @throws Exception
	 */
	public double AngDist(Vector3d v1, Vector3d v2) throws Exception {
		double dist = 0.;
		double aligned = 0.999;
		/* Normalize both vectors */
		Vector3d r1 = new Vector3d(v1.x,v1.y,v1.z);
		Vector3d r2 = new Vector3d(v2.x,v2.y,v2.z);
		r1.normalize();
		r2.normalize();
		double sprod = r1.dot(r2);
		/* This takes care about the bug in vecmath method from java3d project */
		if (sprod > aligned) { // almost aligned
			r1.sub(r2);
			double diff = r1.length();
			dist = 2.0 * Math.asin(diff / 2.0);

		} else if (sprod < -aligned) {
			r1.add(r2);
			double diff = r1.length();
			dist = Math.PI - 2.0 * Math.asin(diff / 2.0);
		} else {
			dist = r1.angle(r2);
		}
		return dist;
	}

	/**
	 * calculates a vector production of two vectors.
	 * 
	 * @param v1 
	 *            Vectror containing 3 elements of Number type
	 * @param v2 
	 *            Vector containing 3 elements of Number type
	 * @return Vector of 3 Objects of Double type
	 * @throws Exception
	 */
	public Vector VectProd(Vector v1, Vector v2) throws Exception {
		Vector res = new Vector();
//
		double[] v1_element = new double[3];
		double[] v2_element = new double[3];
		for (int i = 0; i < 3; i++) {
			if (v1.get(i).getClass().isInstance(Number.class)) {
				v1_element[i] = ((Number) v1.get(i)).doubleValue();
			} else {
				throw new Exception();
			}
			if (v2.get(i).getClass().isInstance(Number.class)) {
				v2_element[i] = ((Number) v2.get(i)).doubleValue();
			} else {
				throw new Exception();
			}

		}

		Double value = new Double(v1_element[1] * v2_element[2] - v1_element[2]
				* v2_element[1]);
		res.add( value);
		value = new Double(v1_element[1] * v2_element[2] - v1_element[2]
				* v2_element[1]);
		res.add( value);
		value = new Double(v1_element[1] * v2_element[2] - v1_element[2]
				* v2_element[1]);
		res.add( value);
		return res;
	}

	/**
	 * calculates a dot product (inner product) of two 3D vectors
	 *  the result is double
	 * 
	 * @param v1 
	 *            3d Vector of Number Objects (Double, long .. )
	 * @param v2 
	 *            3d Vector
	 * @return  double
	 * @throws Exception
	 */
	public double dotProduct(Vector3d v1, Vector3d v2) throws Exception {

		double prod = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;

		return prod;
	}

	/**
	 * calculate cross product of two vectors
	 * 
	 * @param v1 
	 *            Vector3d
	 * @param v2 
	 *            Vector3d
	 * @return  Vector3d result of the product
	 */
	public Vector3d crossProduct(Vector3d v1, Vector3d v2) {
		Vector3d res = new Vector3d(0., 0., 0.);
		double x = v1.y * v2.z - v1.z * v2.y;
		double y = v1.z * v2.x - v1.x * v2.z;
		double z = v1.x * v2.y - v1.y * v2.x;
		res.x = x;
		res.y = y;
		res.z = z;
		return res;
	}
    /**
     * calculates angular resolution of the pixel map
     * in arc seconds.
     * @param nside
     * @return double resolution in arcsec
     */
    public double PixRes(long nside) {
        double res = 0.;
        double degrad = Math.toDegrees(1.0);
        double skyArea = 4.*PI*degrad*degrad; // 4PI steredian in deg^2
        double arcSecArea = skyArea*3600.*3600.;  // 4PI steredian in (arcSec^2)
        long npixels = 12*nside*nside;
        res = arcSecArea/npixels;       // area per pixel
        res = Math.sqrt(res);           // angular size of the pixel arcsec
        return res;
    }
    /**
     * calculate requared nside given pixel size in arcsec
     * @param pixsize in arcsec
     * @return long nside parameter
     */
    public long GetNSide(double pixsize) {
    	long[] nsidelist = { 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048,
				4096, 8192, 16384, 32768, 65536, 131072, 262144, 524288, 1048576 };
    	long res = 0;
    	double pixelArea = pixsize*pixsize;
    	double degrad = Math.toDegrees(1.);
    	double skyArea = 4.*PI*degrad*degrad*3600.*3600.;
    	long npixels = (long) (skyArea/pixelArea);
    	long nsidesq = npixels/12;
    	long nside_req = (long) Math.sqrt(nsidesq);
    	long mindiff = ns_max;
    	int indmin = 0;
    	for (int i=0; i<nsidelist.length; i++) {
    		if (Math.abs(nside_req - nsidelist[i]) <= mindiff) {
    			mindiff = Math.abs(nside_req - nsidelist[i]);
    			res = nsidelist[i];
    			indmin = i;
    		}
    		if ((nside_req > res) && (nside_req < ns_max)) res = nsidelist[indmin+1];
    	   	if (nside_req > ns_max ) {
        		System.out.println("nside cannot be bigger than "+ns_max);
        		return ns_max;
        	}
    		
    	}
    	return res;
    }
    /**
     * returns polar coordinates in radians given ra, dec in degrees
     * @param radec double array containing ra,dec in degrees
     * @return res double array containing theta and phi in radians
     *             res[0] = theta res[1] = phi
     */
    public double[] RaDecToPolar(double[] radec) {
    	double[] res = {0.0,0.0};
    	
			double ra =  radec[0];
			double dec =  radec[1];
			double theta = PI/2. - Math.toRadians(dec);
			double phi = Math.toRadians(ra);
			res[0] = theta;
			res[1] = phi;
    	
    	return res;
    }
    /**
     * returns ra, dec in degrees given polar coordinates in radians
     * @param polar double array polar[0] = phi in radians
     *                           polar[1] = theta in radians
     * @return double array radec radec[0] = ra in degrees
     *                radec[1] = dec in degrees
     */
    public double[] PolarToRaDec(double[] polar) {
    	double[] radec = {0.0,0.0};
			double phi =  polar[1];
			double theta = polar[0];
			double dec = Math.toDegrees(PI/2. - theta);
			double ra = Math.toDegrees(phi);
			radec[0] = ra;
			radec[1] = dec;
    	
    	return radec;
    }
    
    /**
     * returns polar coordinates of a point on unit sphere given Cartesian coordinates
     * @param x - Cartesian coordinate x of a point on unit sphere
     * @param y - y coordinate
     * @param z - z coordinate
     * @return double [] theta,phi
     */
    public double[] xyzToPolar(double x, double y, double z) {
    	double[] res;
    	Vector3d vv = new Vector3d(x,y,z);
    	res = Vect2Ang(vv);
    	return res;
    }

	/**
         * Returns singleton instance.
	 */
	public static PixTools getInstance() {
		return pixTools;
	}
}
