package gov.fnal.eag.healpix.test;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import gov.fnal.eag.healpix.BitManipulation;
import gov.fnal.eag.healpix.PixTools;

import javax.vecmath.Vector3d;

import junit.framework.TestCase;

/**
 * test suit for PixTools using static constructor
 * 
 * @author N. Kuropatkin
 * 
 */

public class PixToolsStaticTest extends TestCase {
    private static PixTools pixTools = new PixTools();

    /**
     *  test MODULO function
     */
    public void testMODULO() {
//
        double A = 8.;
        double B = 5.;
        double res = BitManipulation.MODULO(A, B);
        assertEquals("modulo = " + res, 3., res, 1e-10);
        A = -8.;
        B = 5.;
        res = BitManipulation.MODULO(A, B);
        assertEquals("modulo = " + res, 2., res, 1e-10);
        A = 8.;
        B = -5.;
        res = BitManipulation.MODULO(A, B);
        assertEquals("modulo = " + res, -2., res, 1e-10);
        A = -8.;
        B = -5.;
        res = BitManipulation.MODULO(A, B);
        assertEquals("modulo = " + res, -3., res, 1e-10);
        System.out.println(" test MODULO is done");
    }

    /**
     * tests angular distance calculation
     * @throws Exception
     */
    public void testAngDist() throws Exception {
        Vector3d v1 = new Vector3d(1.5, 1.6, 2.0);
        Vector3d v2 = new Vector3d(1.2, 1.0, 1.5);

        double res1 = pixTools.AngDist(v1, v2);
        double res2 = v1.angle(v2);
        System.out.println("res1 = " + res1 + " res2=" + res2);
        assertEquals("angular Distance=" + res2, 1.0, res1 / res2, 1e-10);
 

        Vector3d v3 = new Vector3d(1.5, 1.6, 0.);
        Vector3d v4 = new Vector3d(1.5, 1.601, 0.);
        double res3 = pixTools.AngDist(v3, v4);
        double res4 = v3.angle(v4);
//
        System.out.println("res3 = " + res3 + " res4=" + res4);
        assertEquals("angular Distance=" + res4, 1.,  res3/ res4, 1e-3);
        Vector3d v5 = new Vector3d(1.5, 1.6, 0.);
        Vector3d v6 = new Vector3d(-1.5, -1.75, 0.);
        double res5 = pixTools.AngDist(v5, v6);
        double res6 = v5.angle(v6);
        System.out.println("res5 = " + res5 + " res6=" + res6);
        assertEquals("angular Distance=" + res6, 1.0, res5 / res6, 1e-10);
        System.out.println(" test of AngDist is done");
    }

    /**
     * @throws Exception
     */
    public void testSurfaceTriangle() throws Exception {
        Vector3d v1 = new Vector3d(1.0, 0.0, 0.0);
        Vector3d v2 = new Vector3d(0.0, 1.0, 0.0);
        Vector3d v3 = new Vector3d(0.0, 0.0, 1.0);
        PixTools pt = new PixTools();
        double res = pixTools.SurfaceTriangle(v1, v2, v3);
        System.out.println("Triangle surface is=" + res / Math.PI
                + " steredians");
        assertEquals("Triangle surface=" + res, 0.5, res / Math.PI, 1e-10);
        System.out.println(" test of SurfaceTriangle is done");
    }

    /**
     * tests calculation of npixels from nside
     */
    public void testNside2Npix() {
        int nside = 1;
        int npix = 0;
        npix = (int) pixTools.Nside2Npix(nside);
        assertEquals("Npix=" + npix, 12, npix, 1e-10);
        nside = 2;
        npix = (int) pixTools.Nside2Npix(nside);
        assertEquals("Npix=" + npix, 48, npix, 1e-10);
    }

    /**
     * tests calculation of nsides from npixels
     */
    public void testNpix2Nside() {
        int npix = 12;
        int nside = 0;
        nside = (int) pixTools.Npix2Nside(npix);

        double pixSize1 = pixTools.PixRes(65536);
        long nside1 = pixTools.GetNSide(pixSize1);
        assertEquals("Nside=" + nside1, 65536, nside1, 1e-10);
        
        assertEquals("Nside=" + nside, 1, nside, 1e-10);

    }

    /**
     * test of directional angles calculation
     */
    public void testVec2Ang() {
        double PI = Math.PI;
        Vector3d v = new Vector3d(0.0, 1.0, 0.0);
        double[] ang_tup = { 0., 0. };
        ang_tup = pixTools.Vect2Ang(v);
        System.out.println(" Theta=" + ang_tup[0] / PI + " Phi=" + ang_tup[1]
                / PI);
        assertEquals("Theta=" + ang_tup[0], 0.5, ang_tup[0] / PI, 1e-10);
        assertEquals("Phi=" + ang_tup[1], 0.5, ang_tup[1] / PI, 1e-10);
        v = new Vector3d(1.0, 0.0, 0.0);
        ang_tup = pixTools.Vect2Ang(v);
        assertEquals("phi=" + ang_tup[1], 0., ang_tup[1] / PI, 1e-10);
        System.out.println(" test Vect2Ang is done");
    }

    /**
     * tests calculation of pixel from polar angles
     * in ring schema of pixalization
     * @throws Exception
     */
    public void testAng2Pix() throws Exception {
        System.out.println(" Test ang2pix_ring ___________________");
        double PI = Math.PI;
        long pix = -1;
        double theta = PI / 2. - 0.2;
        double phi = PI / 2. ; 
        long nside = 4;
        try {
            pix =  pixTools.ang2pix_ring(nside,theta, phi);
        } catch (Exception e) {
            e.printStackTrace();
        }
        Vector3d v = pixTools.Ang2Vec(theta,phi);
        long pix1 =  pixTools.vect2pix_ring(nside,v);
        assertEquals("pix=" + pix, pix1, pix, 1e-10);
        assertEquals("pix=" + pix, 76, pix, 1e-10);

        long pix2 =  pixTools.ang2pix_nest(nside,theta,phi);
        long pix3 =  pixTools.vect2pix_nest(nside,v);
        assertEquals("pix2=" + pix2, pix3, pix2, 1e-10);
        assertEquals("pix2=" + pix2, 92, pix2, 1e-10);


        double[] radec = pixTools.pix2ang_ring(nside,76);
        assertEquals("theta=" + theta, theta, radec[0], 4e-2);
        assertEquals("phi=" + phi, radec[1], phi, 1e-2);
        double[] radec1 = pixTools.pix2ang_nest(nside,92);
        System.out.println("theta="+radec1[0]+" phi="+radec1[1]);
        assertEquals("theta=" + theta, theta, radec1[0], 4e-2);
        assertEquals("phi=" + phi, radec1[1], phi, 1e-2);   
        System.out.println(" test Ang2Pix is done");
    }

    /**
     * tests calculation of unit vector from polar angles
     * @throws Exception
     */
    public void testAng2Vect() throws Exception {
        System.out.println(" Start test Ang2Vect----------------");
        double PI = Math.PI;
        double theta = PI / 2.;
        double phi = PI / 2;
        Vector3d v = pixTools.Ang2Vec(theta, phi);
        System.out.println("Vector x=" + v.x + " y=" + v.y + " z=" + v.z);
        assertEquals("x=" + v.x, 0., v.x, 1e-10);
        assertEquals("y=" + v.y, 1., v.y, 1e-10);
        assertEquals("z=" + v.z, 0., v.z, 1e-10);
        System.out.println(" test Ang2Vect is done");
    }

    /**
     * tests calculation of ring number from z coordinate
     * @throws Exception
     */
    public void testRingNum() throws Exception {
        double z = 0.25;
        int nside = 1;
        System.out.println("Start test RingNum !!!!!!!!!!!!!!!!!!!!");
        int nring = (int) pixTools.RingNum(nside, z);
        System.out.println("z=" + z + " ring number =" + nring);
        assertEquals("z=" + z, 2, nring, 1e-10);
        z = -0.25;
        nring = (int) pixTools.RingNum(nside, z);
        assertEquals("z=" + z, 2, nring, 1e-10);
        z = 0.8;
        nring = (int) pixTools.RingNum(nside, z);
        assertEquals("z=" + z, 1, nring, 1e-10);
        z = -0.8;
        nring = (int) pixTools.RingNum(nside, z);
        assertEquals("z=" + z, 3, nring, 1e-10);
        System.out.println(" test RingNum is done");
        nside = 4;
        int pix = 3;
        Vector3d v = pixTools.pix2vect_ring(nside, pix);
        z = v.z;
        nring = (int) pixTools.RingNum(nside, z);
        assertEquals("z=" + z, 1, nring, 1e-10);
        pix = 11;
        v = pixTools.pix2vect_ring(nside, pix);
        z = v.z;
        nring = (int) pixTools.RingNum(nside, z);
        assertEquals("z=" + z, 2, nring, 1e-10);
        pix = 23;
        v = pixTools.pix2vect_ring(nside, pix);
        z = v.z;
        nring = (int) pixTools.RingNum(nside, z);
        assertEquals("z=" + z, 3, nring, 1e-10);
        pix = 39;
        v = pixTools.pix2vect_ring(nside, pix);
        z = v.z;
        nring = (int) pixTools.RingNum(nside, z);
        assertEquals("z=" + z, 4, nring, 1e-10);
        pix = 55;
        v = pixTools.pix2vect_ring(nside, pix);
        z = v.z;
        nring = (int) pixTools.RingNum(nside, z);
        assertEquals("z=" + z, 5, nring, 1e-10);
        pix = 71;
        v = pixTools.pix2vect_ring(nside, pix);
        z = v.z;
        nring = (int) pixTools.RingNum(nside, z);
        assertEquals("z=" + z, 6, nring, 1e-10);
        pix = 87;
        v = pixTools.pix2vect_ring(nside, pix);
        z = v.z;
        nring = (int) pixTools.RingNum(nside, z);
        assertEquals("z=" + z, 7, nring, 1e-10);
        pix = 103;
        v = pixTools.pix2vect_ring(nside, pix);
        z = v.z;
        nring = (int) pixTools.RingNum(nside, z);
        assertEquals("z=" + z, 8, nring, 1e-10);
        pix = 119;
        v = pixTools.pix2vect_ring(nside, pix);
        z = v.z;
        nring = (int) pixTools.RingNum(nside, z);
        assertEquals("z=" + z, 9, nring, 1e-10);
        pix = 135;
        v = pixTools.pix2vect_ring(nside, pix);
        z = v.z;
        nring = (int) pixTools.RingNum(nside, z);
        assertEquals("z=" + z, 10, nring, 1e-10);
        pix = 151;
        v = pixTools.pix2vect_ring(nside, pix);
        z = v.z;
        nring = (int) pixTools.RingNum(nside, z);
        assertEquals("z=" + z, 11, nring, 1e-10);
        pix = 167;
        v = pixTools.pix2vect_ring(nside, pix);
        z = v.z;
        nring = (int) pixTools.RingNum(nside, z);
        assertEquals("z=" + z, 12, nring, 1e-10);
        pix = 169;
        v = pixTools.pix2vect_ring(nside, pix);
        z = v.z;
        nring = (int) pixTools.RingNum(nside, z);
        assertEquals("z=" + z, 13, nring, 1e-10);
        pix = 180;
        v = pixTools.pix2vect_ring(nside, pix);
        z = v.z;
        nring = (int) pixTools.RingNum(nside, z);
        assertEquals("z=" + z, 14, nring, 1e-10);
        System.out.println("End test RingNum");
    }

    /**
     * tests conversion from nest schema pixel to ring schema pixel
     * @throws Exception
     */
    public void testNest2Ring() throws Exception {
        int ipnest = 3;
        int nside = 2;
        int ipring = (int) pixTools.nest2ring(nside, ipnest);
        assertEquals("ipring=" + ipring, 0, ipring, 1e-10);
        ipnest = 0;
        nside = 2;
        ipring = (int) pixTools.nest2ring(nside, ipnest);
        assertEquals("ipring=" + ipring, 13, ipring, 1e-10);
        ipnest = 18;
        nside = 2;
        ipring = (int) pixTools.nest2ring(nside, ipnest);
        assertEquals("ipring=" + ipring, 27, ipring, 1e-10);
        ipnest = 23;
        nside = 2;
        ipring = (int) pixTools.nest2ring(nside, ipnest);
        assertEquals("ipring=" + ipring, 14, ipring, 1e-10);
        ipnest = 5;
        nside = 4;
        ipring = (int) pixTools.nest2ring(nside, ipnest);
        assertEquals("ipring = " + ipring, 27, ipring, 1e-10);
        System.out.println(" test Nest2Ring is done");
    }

    /**
     * tests conversion from ring schema pixel to nest schema pixel
     * @throws Exception
     */
    public void testRing2Nest() throws Exception {
        System.out.println(" start test Ring2Nest !!!!!!!!!!!!!!!!!!!!!!");
        int ipring = 0;
        int nside = 2;

        int ipnest = (int) pixTools.ring2nest(nside, ipring);
        assertEquals("ipnest=" + ipnest, 3, ipnest, 1e-10);
        ipring = 13;
        nside = 2;
        ipnest = (int) pixTools.ring2nest(nside, ipring);
        assertEquals("ipnest=" + ipnest, 0, ipnest, 1e-10);
        ipring = 27;
        nside = 2;
        ipnest = (int) pixTools.ring2nest(nside, ipring);
        assertEquals("ipnest=" + ipnest, 18, ipnest, 1e-10);
        ipring = 14;
        nside = 2;
        ipnest = (int) pixTools.ring2nest(nside, ipring);
        assertEquals("ipnest=" + ipnest, 23, ipnest, 1e-10);
        ipring = 27;
        nside = 4;
        ipnest = (int) pixTools.ring2nest(nside, ipring);
        assertEquals("ipnest = " + ipnest, 5, ipnest, 1e-10);
        ipring = 83;
        nside = 4;
        ipnest = (int) pixTools.ring2nest(nside, ipring);
        assertEquals("ipnest = " + ipnest, 123, ipnest, 1e-10);
        System.out.println(" test Ring2Nest is done");
    }

    /**
     * tests Next_In_Line method for the nest schema
     * @throws Exception
     */
    public void testNext_In_Line_Nest() throws Exception {
        int ipix = 0;
        int nside = 2;
        int ipnext = (int) pixTools.next_in_line_nest(nside, ipix);
        assertEquals("ipnext=" + ipnext, 23, ipnext, 1e-10);
        ipix = 1;
        nside = 2;
        ipnext = (int) pixTools.next_in_line_nest(nside, ipix);
        assertEquals("ipnext=" + ipnext, 6, ipnext, 1e-10);
        ipix = 4;
        nside = 2;
        ipnext = (int) pixTools.next_in_line_nest(nside, ipix);
        assertEquals("ipnext=" + ipnext, 27, ipnext, 1e-10);
        ipix = 27;
        nside = 2;
        ipnext = (int) pixTools.next_in_line_nest(nside, ipix);
        assertEquals("ipnext=" + ipnext, 8, ipnext, 1e-10);
        ipix = 12;
        nside = 2;
        ipnext = (int) pixTools.next_in_line_nest(nside, ipix);
        assertEquals("ipnext=" + ipnext, 19, ipnext, 1e-10);
        ipix = 118;
        nside = 4;
        ipnext = (int) pixTools.next_in_line_nest(nside, ipix);
        assertEquals("ipnext = " + ipnext, 117, ipnext, 1e-10);
        System.out.println(" test next_in_line_nest is done");
    }

    /**
     * tests InRing method
     * @throws Exception
     */
    public void testInRing() throws Exception {
        System.out.println(" Start test InRing !!!!!!!!!!!!!!!!!!!!!!!!!");
        int[] nestComp = { 19, 0, 23, 4, 27, 8, 31, 12 };
        long[] ringHi = {17, 18, 19, 12, 13 };
        long[] nestHi = { 8, 31,  12, 19, 0 };
        long[] ringLow = {19, 12, 13, 14, 15 };
        long[] nestLow = { 12, 19, 0, 23, 4 };
        double PI = Math.PI;
        int nside = 2;
        boolean nest = false;
        int iz = 3;
        double phi = PI;
        double dphi = PI;
        ArrayList ring = pixTools.InRing(nside, iz, phi, dphi, nest);
        for (int i = 0; i < ring.size(); i++) {
            assertEquals("ipnext = " + ((Long) ring.get(i)).longValue(),
                    i + 12, ((Long) ring.get(i)).longValue(), 1e-10);
        }
        Vector3d v = new Vector3d(1., 0., 0.);
        double[] ang_tup = { 0., 0. };
        ang_tup = pixTools.Vect2Ang(v);
        phi = ang_tup[1]/PI;
        ring = pixTools.InRing(nside, iz, phi, dphi, nest);
        for (int i = 0; i < ring.size(); i++) {
            assertEquals("ipnext = " + ((Long) ring.get(i)).longValue(),
                    i + 12, ((Long) ring.get(i)).longValue(), 1e-10);
        }
        Vector3d v1 = new Vector3d(-1., 0., 0.);
        
        ang_tup = pixTools.Vect2Ang(v1);
        phi = ang_tup[1]/PI;
        ring = pixTools.InRing(nside, iz, phi, dphi, nest);
        for (int i = 0; i < ring.size(); i++) {
            assertEquals("ipnext = " + ((Long) ring.get(i)).longValue(),
                    i + 12, ((Long) ring.get(i)).longValue(), 1e-10);
        }
        phi = 1.75*PI;
        dphi = 0.5*PI;
        ring = pixTools.InRing(nside, iz, phi, dphi, nest);
        for (int i = 0; i < ring.size(); i++) {
            assertEquals("ipnext = " + ((Long) ring.get(i)).longValue(),
                    ringHi[i], ((Long) ring.get(i)).longValue(), 1e-10);

        }

        phi = 1.75*PI;
        dphi = 0.5*PI;
        nest = true;
        ring = pixTools.InRing(nside, iz, phi, dphi, nest);
        for (int i = 0; i < ring.size(); i++) {
            assertEquals("ipnext = " + ((Long) ring.get(i)).longValue(),
                    nestHi[i], ((Long) ring.get(i)).longValue(), 1e-10);

        }   
        phi = 0.25*PI;
        dphi = 0.5*PI;
        nest = false;
        ring = pixTools.InRing(nside, iz, phi, dphi, nest);
        for (int i = 0; i < ring.size(); i++) {
            assertEquals("ipnext = " + ((Long) ring.get(i)).longValue(),
                    ringLow[i], ((Long) ring.get(i)).longValue(), 1e-10);

        }           
        phi = 0.25*PI;
        dphi = 0.5*PI;
        nest = true;
        ring = pixTools.InRing(nside, iz, phi, dphi, nest);
        for (int i = 0; i < ring.size(); i++) {
            assertEquals("ipnext = " + ((Long) ring.get(i)).longValue(),
                    nestLow[i], ((Long) ring.get(i)).longValue(), 1e-10);

        }   
        
        nest = true;
        dphi = PI;
        ring = pixTools.InRing(nside, iz, phi, dphi, nest);
        for (int i = 0; i < ring.size(); i++) {
            assertEquals("ipnext = " + ((Long) ring.get(i)).longValue(),
                    nestComp[i], ((Long) ring.get(i)).longValue(), 1e-10);
        }
        nest = false;
        nside = 4;
        phi = 2.1598449493429825;
        iz = 8;
        dphi = 0.5890486225480867;
        //      System.out.println(" iz="+iz+" phi="+phi+" dphi="+dphi);
        ring = pixTools.InRing(nside, iz, phi, dphi, nest);
        //      for (int i = 0; i<ring.size(); i++) {
        //          System.out.println("ipnext = "+ ((Integer)ring.get(i)).intValue());
        //      }
        nest = false;
        nside = 4;
        dphi = 0. * PI;
        iz = 8;
        phi = 2.1598449493429825;
        //      System.out.println(" iz="+iz+" phi="+phi+" dphi="+dphi);
        ring = pixTools.InRing(nside, iz, phi, dphi, nest);
        //      for (int i = 0; i<ring.size(); i++) {
        //          System.out.println("ipnext = "+ ((Integer)ring.get(i)).intValue());
        //      }
        System.out.println(" test InRing is done");
    }

    /**
     * tests Neighbour's method for nest schema of the pixelization 
     * @throws Exception
     */
    public void testNeighbours_Nest() throws Exception {
        System.out.println(" Start test Neighbours_Nest !!!!!!!!!!!!!!!!!");
        long nside = 2;
        long ipix = 25;
        long[] test17 = { 34, 16, 18, 19, 2, 0, 22, 35 };
        long[] test32 = { 40, 44, 45, 34, 35, 33, 38, 36 };
        long[] test3 = { 0, 2, 13, 15, 11, 7, 6, 1 };
        long[] test25 = { 42, 24, 26, 27, 10, 8, 30, 43 };
        //
        ArrayList npixList = new ArrayList();
        npixList = pixTools.neighbours_nest(nside, ipix);
        for (int i = 0; i < npixList.size(); i++) {
            assertEquals("ip = " + ((Long) npixList.get( i)).longValue(),
                    test25[ i], ((Long) npixList.get(  i)).longValue(), 1e-10);
        }
        ipix = 17;

        npixList = pixTools.neighbours_nest(nside, ipix);
        for (int i = 0; i < npixList.size(); i++) {
            assertEquals("ip = " + ((Long) npixList.get(i)).longValue(),
                    test17[i], ((Long) npixList.get(i)).longValue(), 1e-10);
        }
        ipix = 32;

        npixList = pixTools.neighbours_nest(nside, ipix);
        for (int i = 0; i < npixList.size(); i++) {
            assertEquals("ip = " + ((Long) npixList.get(i)).longValue(),
                    test32[i], ((Long) npixList.get(i)).longValue(), 1e-10);
        }
        ipix = 3;

        npixList = pixTools.neighbours_nest(nside, ipix);
        for (int i = 0; i < npixList.size(); i++) {
            assertEquals("ip = " + ((Long) npixList.get(i)).longValue(),
                    test3[i], ((Long) npixList.get(i)).longValue(), 1e-10);
        }
        System.out.println(" test NeighboursNest is done");
    }

    /**
     * tests intrs_intrv method
     * @throws Exception
     */
    public void testIntrs_Intrv() throws Exception {
        System.out.println(" test intrs_intrv !!!!!!!!!!!!!!!!!!!!!!!!!!!");
        double[] d1 = { 1.0, 9.0 };
        double[] d2 = { 3.0, 16.0 };
        double[] di;
        //      System.out.println("Case "+d1[0]+" "+d1[1]+" | "+d2[0]+" "+d2[1]);
        di = pixTools.intrs_intrv(d1, d2);
        //      System.out.println("Result "+di[0]+" - "+di[1]);
        int n12 = di.length / 2;
        assertEquals("n12 = " + n12, 1, n12, 1e-6);
        assertEquals("di[0] = " + di[0], 3.0, di[0], 1e-6);
        assertEquals("di[1] = " + di[1], 9.0, di[1], 1e-6);
        d1 = new double[] { 0.537, 4.356 };
        d2 = new double[] { 3.356, 0.8 };
        //      System.out.println("Case "+d1[0]+" "+d1[1]+" | "+d2[0]+" "+d2[1]);
        di = pixTools.intrs_intrv(d1, d2);
        n12 = di.length / 2;
        assertEquals("n12 = " + n12, 2, n12, 1e-6);
        assertEquals("di[0] = " + di[0], 0.537, di[0], 1e-6);
        assertEquals("di[1] = " + di[1], 0.8, di[1], 1e-6);
        assertEquals("di[2] = " + di[2], 3.356, di[2], 1e-6);
        assertEquals("di[1] = " + di[3], 4.356, di[3], 1e-6);

        d1 = new double[] { 2.356194490092345, 2.356194490292345 };
        d2 = new double[] { 1.251567, 4.17 };
        //      System.out.println("Case "+d1[0]+" "+d1[1]+" | "+d2[0]+" "+d2[1]);
        di = pixTools.intrs_intrv(d1, d2);
        n12 = di.length / 2;
        assertEquals("n12 = " + n12, 1, n12, 1e-6);
        assertEquals("di[0] = " + di[0], 2.35619449009, di[0], 1e-6);
        assertEquals("di[1] = " + di[1], 2.35619449029, di[1], 1e-6);

        System.out.println(" test intrs_intrv is done");
    }

    /**
     * tests conversion from pixel number to vector
     * @throws Exception
     */
    public void testPix2Vect_ring() throws Exception {
        System.out.println("Start test Pix2Vect_ring !!!!!!!!!!!!!!!!!!!");
        double TWOPI = 2.0 * Math.PI;
        int nside = 2;
        int ipix = 0;
        Vector3d v1 = new Vector3d(0., 0., 0.);
        v1 = pixTools.pix2vect_ring(nside, ipix);
        assertEquals("v1.z = " + v1.z, 1.0, v1.z, 1e-1);

        ipix = 20;
        Vector3d v2 = new Vector3d(0., 0., 0.);
        v2 = pixTools.pix2vect_ring(nside, ipix);
        assertEquals("v2.x = " + v2.x, 1.0, v2.x, 1e-1);
        assertEquals("v2.z = " + v2.z, 0.0, v2.z, 1e-1);
        ipix = 22;
        Vector3d v3 = new Vector3d();
        v3 = pixTools.pix2vect_ring(nside, ipix);
        assertEquals("v3.y = " + v3.y, 1.0, v3.y, 1e-1);
        assertEquals("v3.z = " + v3.z, 0.0, v3.z, 1e-1);
        //      System.out.println("Vector3 x="+v3.x+" y="+v3.y+" z="+v3.z);
        ipix = 95;
        nside = 4;
        v1 = pixTools.pix2vect_ring(nside, ipix);
        v1.normalize();
        double phi1 = Math.atan2(v1.y, v1.x);
        double[] tetphi = new double[2];
        tetphi = pixTools.pix2ang_ring(nside, ipix);
        assertEquals("phi = " + phi1, 0.0, Math.abs(phi1 - tetphi[1]), 1e-10);
        ipix = 26;
        nside = 4;
        v1 = pixTools.pix2vect_ring(nside, ipix);
        v1.normalize();
        phi1 = Math.atan2(v1.y, v1.x);
        if (phi1 < 0.)
            phi1 += TWOPI;
        tetphi = new double[2];
        tetphi = pixTools.pix2ang_ring(nside, ipix);
        assertEquals("phi = " + phi1, 0.0, Math.abs(phi1 - tetphi[1]), 1e-10);
        System.out.println("------------------------------------------");
        System.out.println(" test pix2vect_ring is done");
    }

    /**
     * tests conversion from pixel number to vector
     * @throws Exception
     */
    public void testPix2Vect_nest() throws Exception {
        double TWOPI = 2.0 * Math.PI;
        int nside = 2;
        int ipix = 3;
        System.out.println(" Start test Pix2Vect_nest !!!!!!!!!!!!!!");
        Vector3d v1 = new Vector3d(0., 0., 0.);
        v1 = pixTools.pix2vect_nest(nside, ipix);
        assertEquals("v1.z = " + v1.z, 1.0, v1.z, 1e-1);
        ipix = 17;
        Vector3d v2 = new Vector3d(0., 0., 0.);
        v2 = pixTools.pix2vect_nest(nside, ipix);
        assertEquals("v2.x = " + v2.x, 1.0, v2.x, 1e-1);
        assertEquals("v2.z = " + v2.z, 0.0, v2.z, 1e-1);

        ipix = 21;
        Vector3d v3 = new Vector3d();
        v3 = pixTools.pix2vect_nest(nside, ipix);
        assertEquals("v3.y = " + v3.y, 1.0, v3.y, 1e-1);
        assertEquals("v3.z = " + v3.z, 0.0, v3.z, 1e-1);
        nside = 4;
        ipix = 105;
        v1 = pixTools.pix2vect_nest(nside, ipix);
        v1.normalize();
        double phi1 = Math.atan2(v1.y, v1.x);
        if (phi1 < 0.)
            phi1 += TWOPI;
        double[] tetphi = new double[2];
        tetphi = pixTools.pix2ang_nest(nside, ipix);
        assertEquals("phi = " + phi1, 0.0, Math.abs(phi1 - tetphi[1]), 1e-10);
        nside = 4;
        ipix = 74;
        v1 = pixTools.pix2vect_nest(nside, ipix);
        v1.normalize();
        phi1 = Math.atan2(v1.y, v1.x);
        if (phi1 < 0.)
            phi1 += TWOPI;
        tetphi = new double[2];
        tetphi = pixTools.pix2ang_nest(nside, ipix);
        assertEquals("phi = " + phi1, 0.0, Math.abs(phi1 - tetphi[1]), 1e-10);

        System.out.println(" test pix2vect_nest is done");
        System.out.println("-------------------------------------------");
    }

    /**
     * tests conversion from vector to pixel number
     * @throws Exception
     */
    public void testVect2Pix_ring() throws Exception {
        System.out.println("Start test Vect2Pix_ring !!!!!!!!!!!!!!!!!!!");

        long nside = 4;
        long ipix = 83;
        long respix = 0;
        Vector3d v1 = new Vector3d(0., 0., 0.);
        v1 = pixTools.pix2vect_ring(nside, ipix);
        respix = (int) pixTools.vect2pix_ring(nside, v1);
        assertEquals("respix = " + respix, 83, respix, 1e-10);
        // Hi resolution test
        long nside1 = 1 << 20;
        long maxpix= pixTools.Nside2Npix(nside1);
        System.out.println("nside="+nside1+" maxpix="+maxpix);
        Vector3d v2 = new Vector3d( -0.704, 0.580, 0.408);
        respix = pixTools.vect2pix_ring(nside1, v2);  // ring pixel
        long respixN = pixTools.ring2nest(nside1, respix); // convert to nest
        long respixNC = pixTools.vect2pix_nest(nside1,v2); // nest pixel from the same vector
        long respixR = pixTools.nest2ring(nside1, respixN); // convert pixel 
        System.out.println(" orig="+respix+" doubleT="+respixR+" nest="+respixN+" correct nest="+respixNC);
        assertEquals("ringpix = " + respix, respix, respixR);
        assertEquals("nestpix = " + respixNC, respixNC, respixN);
        System.out.println("------------------------------------------");
        System.out.println(" test vect2pix_ring is done");
    }

    /**
     * tests conversion from vector to pixel number
     * @throws Exception
     */
    public void testVect2Pix_nest() throws Exception {
        System.out.println("Start test Vect2Pix_nest !!!!!!!!!!!!!!!!!!!");
        long nside = 4;
        long ipix = 83;
        long respix = 0;
        Vector3d v1 = new Vector3d(0., 0., 0.);
        v1 = pixTools.pix2vect_ring(nside, ipix);
        respix = (int) pixTools.vect2pix_nest(nside, v1);
        assertEquals("respix = " + respix, 123, respix, 1e-10);
        //
        long nside1 = 1 << 20;
        long maxpix= pixTools.Nside2Npix(nside1);
        System.out.println("nside="+nside1+" maxpix="+maxpix);
        Vector3d v2 = new Vector3d( -0.704, 0.580, 0.408);
        respix = pixTools.vect2pix_nest(nside1, v2);
        long respixRC = pixTools.vect2pix_ring(nside1, v2);
        long respixR = pixTools.nest2ring(nside1, respix);
        long respixN = pixTools.ring2nest(nside1, respixRC);
        System.out.println(" orig="+respix+" doubleT="+respixN+" ring="+respixR+" correct ring="+respixRC);
        assertEquals("ringpix = " + respixRC, respixRC, respixR);
        assertEquals("nestpix = " + respix, respix, respixN);
        System.out.println("------------------------------------------");
        System.out.println(" test vect2pix_nest is done");
    }

    /**
     * tests Query_Strip method
     * @throws Exception
     */
    public void testQuery_Strip() throws Exception {
        System.out.println(" Start test query Strip !!!!!!!!!!!!!!!!");
        int[] pixel1 = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
                16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
                32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
                48, 49, 50, 51, 52, 53, 54, 55 };
        int nside = 4;
        int nest = 0;
        double theta1 = 0.0;
        double theta2 = Math.PI / 4.0 + 0.2;
        ArrayList pixlist;
        pixlist = pixTools.query_strip(nside, theta1, theta2, nest);
        int nlist = pixlist.size();
        for (int i = 0; i < nlist; i++) {
            long ipix = ((Long) pixlist.get(i)).longValue();
            assertEquals("pixel = " + ipix, pixel1[i], ipix, 1e-10);
        }

        System.out.println(" test query_strip is done");

    }

    /**
     * tests Query_Disc method
     * @throws Exception
     */
    public void testQuery_Disc() throws Exception {
        System.out.println(" Start test query_disc !!!!!!!!!!!!!!!!!!!!!");
        long nside = 4;
        int nest = 0;
        long ipix = 0;


          int[] pixel1 = { 45, 46, 60, 61, 62, 77, 78,  92, 93, 94,
                     109, 110,  124, 125, 126, 141, 142
                     };
          


          int[] pixel2 = { 24, 19, 93, 18, 17,  87, 16,86, 85,
                    106,  84, 159,  81, 158, 157, 155, 156
                    };
          

          int[] pixel3 = {52, 79, 49, 78, 77,  75, 76,  74, 73, 70,
                  72, 67,  189, 66, 65, 183, 64
                  };
        int inclusive = 1;
        double radius = Math.PI / 8.0;
        Vector3d v = pixTools.pix2vect_ring(nside, 93);
        ArrayList pixlist;
        pixlist = pixTools.query_disc(nside, v, radius, nest, inclusive);

        int nlist = pixlist.size();

        for (int i = 0; i < nlist; i++) {
            ipix = ((Long) pixlist.get(i)).longValue();
//          System.out.println("i="+i+"pixel="+ipix+" pixel1="+pixel1[i]);
            assertEquals("pixel = " + ipix, pixel1[i], ipix, 1e-10);
        }
        radius = Math.PI/2.0; 
        nest = 0;
        v = new Vector3d(1.0,0.0,0.0);
        ArrayList pixlistT = pixTools.query_disc(2, v, radius, nest, inclusive);

        radius = Math.PI / 8.0;
        v = pixTools.pix2vect_ring(nside, 93);
        nest = 1;
        ArrayList pixlist2 = pixTools.query_disc(nside, v, radius, nest, inclusive);
        
        int nlist2 = pixlist2.size();
        assertEquals("npix="+nlist, nlist, nlist2, 1e-10);
        for (int i = 0; i < nlist2; i++) {
            ipix = ((Long) pixlist2.get(i)).longValue();

            assertEquals("pixel = " + ipix, pixel2[i], ipix, 1e-10);

        }

        v = pixTools.pix2vect_ring(nside, 103);
        nest = 1;
        inclusive = 1;
        ArrayList pixlist3 = pixTools.query_disc(nside, v, radius, nest, inclusive);
        nlist = pixlist3.size();
        assertEquals("npix="+nlist, nlist, pixel3.length, 1e-10);
        
        for (int i = 0; i < pixlist3.size(); i++) {
            ipix = ((Long) pixlist3.get(i)).longValue();

            assertEquals("pixel = " + ipix, pixel3[i], ipix, 1e-10);
        }


        for (int i=0; i<pixel1.length; i++) {
            long ipixR = pixel1[i]; 
            long ipixT = pixTools.ring2nest(nside, ipixR);
            assertEquals("pixel="+ipixT, ipixT, pixel2[i], 1e-10);
            long ipixN = pixTools.nest2ring(nside, ipixT);
            assertEquals("pixel="+ipixN, ipixN, pixel1[i], 1e-10);
        }
        System.out.println(" End query_disk test ______________________________");
    }

    /**
     * tests Query_Triangle method
     * @throws Exception
     */
    public void testQuery_Triangle() throws Exception {
        long nside = 4;
        int nest = 0;
        long ipix = 0;
        int[] pixel1 = { 57, 58, 59, 60, 61, 62, 74, 75, 76, 77, 78, 90, 91,
                92, 93, 107, 108, 109, 123, 124, 140 };
        int[] pixel2 = { 88, 89, 90, 91,  105, 106, 107, 108, 121, 122, 123,
                138, 139, 154 };
        int[] pixel3 = { 49, 64,  80, 81, 95, 96,  112, 113, 
                127, 128, 142, 143, 144, 145 };
        int[] pixel4 = { 36, 52, 53, 67, 68, 69, 83, 84, 85, 86, 98, 99, 100,
                101, 102, 114, 115, 116, 117, 118, 119, 129, 130, 131, 132,
                133, 134, 135 };
        int[] pixel5 = { 58, 127, 56, 126, 125, 50, 123, 124, 119, 48, 122,
                121, 118, 117, 74, 175, 120, 115, 116, 191, 72, 174, 173, 114,
                113, 190, 189, 66 };
        int[] pixel6 = { 110, 123, 124, 125, 140, 141, 156 };
        int[] pixel7 = { 53, 68, 69 };
        long pix1 = 62;
        long pix2 = 57;
        long pix3 = 140;
        System.out.println("Start test Query Triangle !!!!!!!!!!!!!!!!!!!!");
        Vector3d v11 = pixTools.pix2vect_ring(nside, pix1);

        Vector3d v22 = pixTools.pix2vect_ring(nside, pix2);

        Vector3d v33 = pixTools.pix2vect_ring(nside, pix3);

        //      System.out.println("nside="+nside+" triangle pixels "+pix1+" "+pix2+"
        // "+pix3);
        int inclusive = 0;

        ArrayList pixlist;
        pixlist = pixTools.query_triangle(nside, v11, v22, v33, nest, inclusive);

        int nlist = pixlist.size();

        for (int i = 0; i < nlist; i++) {
            ipix = ((Long) pixlist.get(i)).longValue();

            assertEquals("pixel = " + ipix, pixel1[i], ipix, 1e-10);
        }
        pix1 = 92;
        pix2 = 88;
        pix3 = 154;
        v11 = pixTools.pix2vect_ring(nside, pix1);
        v22 = pixTools.pix2vect_ring(nside, pix2);
        v33 = pixTools.pix2vect_ring(nside, pix3);

        inclusive = 0;
        ArrayList pixlist1;
        pixlist1 = pixTools.query_triangle(nside, v11, v22, v33, nest, inclusive);

        nlist = pixlist1.size();

        for (int i = 0; i < nlist; i++) {
            ipix = ((Long) pixlist1.get(i)).longValue();
//          System.out.println(ipix);
            assertEquals("pixel = " + ipix, pixel2[i], ipix, 1e-10);
        }
        pix1 = 49;
        pix2 = 142;
        pix3 = 145;
        v11 = pixTools.pix2vect_ring(nside, pix1);
        v22 = pixTools.pix2vect_ring(nside, pix2);
        v33 = pixTools.pix2vect_ring(nside, pix3);

        inclusive = 0;
        ArrayList pixlist2;
        pixlist2 = pixTools.query_triangle(nside, v11, v22, v33, nest, inclusive);

        nlist = pixlist2.size();

        for (int i = 0; i < nlist; i++) {
            ipix = ((Long) pixlist2.get(i)).longValue();
//          System.out.println(ipix);
            assertEquals("pixel = " + ipix, pixel3[i], ipix, 1e-10);
        }
        pix1 = 36;
        pix2 = 129;
        pix3 = 135;
        v11 = pixTools.pix2vect_ring(nside, pix1);
        v22 = pixTools.pix2vect_ring(nside, pix2);
        v33 = pixTools.pix2vect_ring(nside, pix3);

        inclusive = 0;

        pixlist2 = pixTools.query_triangle(nside, v11, v22, v33, nest, inclusive);

        nlist = pixlist2.size();

        for (int i = 0; i < nlist; i++) {
            ipix = ((Long) pixlist2.get(i)).longValue();
            assertEquals("pixel = " + ipix, pixel4[i], ipix, 1e-10);
        }
        pix1 = 36;
        pix2 = 129;
        pix3 = 135;
        nest = 1;
        v11 = pixTools.pix2vect_ring(nside, pix1);
        v22 = pixTools.pix2vect_ring(nside, pix2);
        v33 = pixTools.pix2vect_ring(nside, pix3);
        inclusive = 0;

        pixlist2 = pixTools.query_triangle(nside, v11, v22, v33, nest, inclusive);

        nlist = pixlist2.size();

        for (int i = 0; i < nlist; i++) {
            ipix = ((Long) pixlist2.get(i)).longValue();
//            System.out.println("ipix="+ipix);
            assertEquals("pixel = " + ipix, pixel5[i], ipix, 1e-10);
        }
        pix1 = 123;
        pix2 = 156;
        pix3 = 110;
        nest = 0;
        v11 = pixTools.pix2vect_ring(nside, pix1);
        v22 = pixTools.pix2vect_ring(nside, pix2);
        v33 = pixTools.pix2vect_ring(nside, pix3);
        inclusive = 0;

        pixlist2 = pixTools.query_triangle(nside, v11, v22, v33, nest, inclusive);

        nlist = pixlist2.size();

        for (int i = 0; i < nlist; i++) {
            ipix = ((Long) pixlist2.get(i)).longValue();
            assertEquals("pixel = " + ipix, pixel6[i], ipix, 1e-10);
            //          System.out.println("i="+i+" pixel#="+ipix);
        }
        pix1 = 69;
        pix2 = 53;
        pix3 = 68;
        nest = 0;
        v11 = pixTools.pix2vect_ring(nside, pix1);
        v22 = pixTools.pix2vect_ring(nside, pix2);
        v33 = pixTools.pix2vect_ring(nside, pix3);
        inclusive = 0;

        pixlist2 = pixTools.query_triangle(nside, v11, v22, v33, nest, inclusive);

        nlist = pixlist2.size();

        for (int i = 0; i < nlist; i++) {
            ipix = ((Long) pixlist2.get(i)).longValue();
            assertEquals("pixel = " + ipix, pixel7[i], ipix, 1e-10);
            //          System.out.println("i="+i+" pixel#="+ipix);
        }
        System.out.println(" test query_triangle is done");

    }

    /**
     * tests Query_Poligon method
     * @throws Exception
     */
    public void testQuery_Polygon() throws Exception {
        long nside = 4;
        int nest = 0;
        long ipix = 0;
        int inclusive = 0;
        int[] result = { 51, 52, 53, 66, 67, 68, 69, 82, 83, 84, 85, 86, 98,
                99, 100, 101, 115, 116, 117 };
        int[] result1 = { 55, 70, 71, 87 };
        int[] result2 = { 137, 152, 153, 168 };
        int[] result3 = { 27, 43, 44, 58, 59, 60, 74, 75, 76, 77, 89, 90, 91,
                92, 93, 105, 106, 107, 108, 109, 110, 121, 122, 123, 124, 125,
                138, 139, 140, 141, 154, 156 };


        System.out.println("Start test query_polygon !!!!!!!!!!!!!!!!!!!!!!");
        ArrayList vlist = new ArrayList();
        Vector3d v = pixTools.pix2vect_ring(nside, 53);
        vlist.add( v);
        v = pixTools.pix2vect_ring(nside, 51);
        vlist.add( v);
        v = pixTools.pix2vect_ring(nside, 82);
        vlist.add( v);
        v = pixTools.pix2vect_ring(nside, 115);
        vlist.add( v);
        v = pixTools.pix2vect_ring(nside, 117);
        vlist.add( v);
        v = pixTools.pix2vect_ring(nside, 86);
        vlist.add( v);

        ArrayList pixlist;
        pixlist = pixTools.query_polygon(nside, vlist, nest, inclusive);
        //      System.out.println(" List size="+pixlist.size());
        int nlist = pixlist.size();
        //      System.out.println(" Pixel list:");
        for (int i = 0; i < nlist; i++) {
            ipix = ((Long) pixlist.get(i)).longValue();
            assertEquals("pixel = " + ipix, result[i], ipix, 1e-10);
            //          System.out.println("i="+i+" pixel # "+ipix);
        }

        /* Yet another test */

        ArrayList vlist1 = new ArrayList();
        v = pixTools.pix2vect_ring(nside, 71);
        vlist1.add( v);
        v = pixTools.pix2vect_ring(nside, 55);
        vlist1.add( v);
        v = pixTools.pix2vect_ring(nside, 70);
        vlist1.add( v);
        v = pixTools.pix2vect_ring(nside, 87);
        vlist1.add( v);
        pixlist = pixTools.query_polygon(nside, vlist1, nest, inclusive);

        nlist = pixlist.size();

        for (int i = 0; i < nlist; i++) {
            ipix = ((Long) pixlist.get(i)).longValue();
            //          System.out.println("i="+i+" pixel # "+ipix);
            assertEquals("pixel = " + ipix, result1[i], ipix, 1e-10);
        }

        /* Yet another test */
        ArrayList vlist2 = new ArrayList();
        v = pixTools.pix2vect_ring(nside, 153);
        vlist2.add( v);
        v = pixTools.pix2vect_ring(nside, 137);
        vlist2.add( v);
        v = pixTools.pix2vect_ring(nside, 152);
        vlist2.add( v);
        v = pixTools.pix2vect_ring(nside, 168);
        vlist2.add( v);
        pixlist = pixTools.query_polygon(nside, vlist2, nest, inclusive);

        nlist = pixlist.size();

        for (int i = 0; i < nlist; i++) {
            ipix = ((Long) pixlist.get(i)).longValue();
            assertEquals("pixel = " + ipix, result2[i], ipix, 1e-10);
            //          System.out.println("i="+i+" pixel # "+ipix);
        }
        /* Yet another test */

        ArrayList vlist3 = new ArrayList();
        v = pixTools.pix2vect_ring(nside, 110);
        vlist3.add( v);
        v = pixTools.pix2vect_ring(nside, 27);
        vlist3.add( v);
        v = pixTools.pix2vect_ring(nside, 105);
        vlist3.add( v);
        v = pixTools.pix2vect_ring(nside, 154);
        vlist3.add( v);
        v = pixTools.pix2vect_ring(nside, 123);
        vlist3.add( v);
        v = pixTools.pix2vect_ring(nside, 156);
        vlist3.add( v);
        pixlist = pixTools.query_polygon(nside, vlist3, nest, inclusive);

        nlist = pixlist.size();

        for (int i = 0; i < nlist; i++) {
            ipix = ((Long) pixlist.get(i)).longValue();
            assertEquals("pixel = " + ipix, result3[i], ipix, 1e-10);
            //                      System.out.println("i="+i+" pixel # "+ipix);
        }
        System.out.println(" test query_polygon is done");

    }
    /**
     * tests MaxResolution method
     */
    public void testMaxResolution() {
        System.out.println(" Start test MaxRes !!!!!!!!!!!!!!!!!!!!!");

        long nside = 1048576;
        double res = pixTools.PixRes(nside);
        System.out.println("Minimum size of the pixel side is "+res+" arcsec.");
        assertEquals("res = " + res, 0.2, res, 1e-1);
        long nsideR = pixTools.GetNSide(res);
        assertEquals("nside = " + nside, nside, nsideR, 1e-1);
        System.out.println(" End of MaxRes test _______________________");
    }
    /**
     * tests QueryDiscResolution method
     * @throws Exception
     */
    public void testQueryDiscRes() throws Exception {
        System.out.println(" Start test DiscRes !!!!!!!!!!!!!!!!!!!!!");

        int nest = 0;
        long ipix = 0;
 
        int inclusive = 0;
        double theta= Math.PI;
        double phi = Math.PI;
        double radius = Math.toRadians(0.2/3600.); //  One arcse
        long nside = pixTools.GetNSide(radius);
        System.out.println(" calculated nside="+nside);
        long cpix = pixTools.ang2pix_ring(nside,theta,phi);
        Vector3d vc = pixTools.pix2vect_ring(nside, cpix);
        ArrayList pixlist;
        pixlist = pixTools.query_disc(nside, vc, radius, nest, inclusive);

        int nlist = pixlist.size();
        for (int i = 0; i < nlist; i++) {
            ipix = ((Long) pixlist.get(i)).longValue();
            Vector3d v = pixTools.pix2vect_ring(nside,ipix);
            double dist = pixTools.AngDist(v,vc);
            assertTrue(dist<=2.*radius);
        }
        cpix = pixTools.ang2pix_nest(nside,theta,phi);
        Vector3d vc1 = pixTools.pix2vect_nest(nside, cpix);
        ArrayList pixlist1;
        nest = 1;
        radius *=4;
        pixlist1 = pixTools.query_disc(nside, vc1, radius, nest, inclusive);
        int nlist1 = pixlist1.size();
        for (int i = 0; i < nlist1; i++) {
            ipix = ((Long) pixlist1.get(i)).longValue();
            Vector3d v = pixTools.pix2vect_nest(nside,ipix);
            double dist = pixTools.AngDist(v,vc1);
            assertTrue(dist<=2.*radius);
        }
        System.out.println(" test query disk  is done -------------------"); 
    }
    /**
     * test Query_disk  check for consistency in the query for RING/NESTED
     */
    public void testQuery_disk2() {
        System.out.println(" Start test query_disk HiRes!!!!!!!!!!!!!!!!!!!!!!!!");
        long nside = 1 << 20 ;
        double res = pixTools.PixRes(nside);
        System.out.println("nside="+nside+" sresolution="+res);
        double radius = Math.toRadians(res/3600.)/2.;
        System.out.println("radius="+radius);
        Vector3d v1 = new Vector3d(-0.704, 0.580, 0.408);
        System.out.println("!!!!!!!!!!!!! NESTED !!!!!!!!!!!");
        ArrayList diskQ = pixTools.query_disc(nside,
                v1,
                radius, 1,1);  // inclusive query at vector point
        assertEquals("npixels = " + diskQ.size(), 8, diskQ.size() , 1e-1);
        long pix1 = pixTools.vect2pix_nest(nside, v1);
        Vector3d v2 = pixTools.pix2vect_nest(nside, pix1);  // vector to pix center
        //
        ArrayList diskQ2 = pixTools.query_disc(nside,
                v2,
                radius, 1,1);  // inclusive with point at pixel center
        assertEquals("npixels = " + diskQ2.size(), 9, diskQ2.size() , 1e-1);

        //
        ArrayList diskQ3 = pixTools.query_disc(nside,
                v2,
                radius, 1,0);  // exclusive with point at pixel center
        assertEquals("npixels = " + diskQ3.size(), 1, diskQ3.size() , 1e-1);

 //   RING schema   
        System.out.println("!!!!!!!!!!!!! RING !!!!!!!!!!!");
        ArrayList diskQ4 = pixTools.query_disc(nside,
                v1,
                radius, 0,1);   // inclusiv at vector point 
        assertEquals("npixels = " + diskQ4.size(), 8, diskQ4.size() , 1e-1);
        //
  
        ArrayList diskQ5 = pixTools.query_disc(nside,
                v2,
                radius, 0, 1);  // inclusive at pixel center
        assertEquals("npixels = " + diskQ5.size(), 9, diskQ5.size() , 1e-1);

//      System.out.println("n pixels in disk5 ="+diskQ5.size());
        ArrayList diskQ6 = pixTools.query_disc(nside,
                v2,
                radius, 0,0);  // exclusive at pixel center
        assertEquals("npixels = " + diskQ6.size(), 1, diskQ6.size() , 1e-1);
//
//  test HiRes conversions
//      
        Vector3d pos = new Vector3d( -0.704, 0.580, 0.408 );

        nside = 1 << 20;
        System.out.println("HiRes transformation tests: nside="+nside);
        ArrayList nestPixels = pixTools.query_disc(nside, pos, radius, 1, 1);
        ArrayList ringPixels = pixTools.query_disc(nside, pos, radius, 0, 1);
        assertEquals(nestPixels.size(), ringPixels.size());
        for(int i=0; i< ringPixels.size(); i++) {
            long iring = ((Number)ringPixels.get(i)).longValue();
            Vector3d cv = pixTools.pix2vect_ring(nside, iring);
            long inest = pixTools.ring2nest(nside, iring);
            long inestC = ((Number)nestPixels.get(i)).longValue();
            Vector3d cvN = pixTools.pix2vect_nest(nside, inestC);
            long iringT = pixTools.nest2ring(nside, inestC);
            assertEquals(iring,iringT);
            assertEquals(inest,inestC);
            assertEquals(" Xv="+cv.x,cv.x,cvN.x,1.e-10);
            assertEquals(" Yv="+cv.y,cv.y,cvN.y,1.e-10);
            assertEquals(" Zv="+cv.z,cv.z,cvN.z,1.e-10);
//          System.out.println(" inest orig="+inestC+" transformed="+inest+" iring orig="+iring+" transf="+iringT);
//          System.out.println("Vector cv vs cvN x="+cv.x+" cvN.x="+cvN.x);
//          System.out.println("Vector cv vs cvN y="+cv.y+" cvN.y="+cvN.y);
//          System.out.println("Vector cv vs cvN z="+cv.z+" cvN.z="+cvN.z);
            double[] tetphiR = pixTools.pix2ang_ring(nside, iring);
            double[] tetphiN= pixTools.pix2ang_nest(nside, inestC);
            assertEquals(" theta="+tetphiR[0],tetphiR[0],tetphiN[0],1.e-10);
            assertEquals(" phi="+tetphiR[1],tetphiR[1],tetphiN[1],1.e-10);
//          System.out.println("theta R vs N "+tetphiR[0]+" "+tetphiN[0]);
//          System.out.println("phi R vs N "+tetphiR[1]+" "+tetphiN[1]);
        }
        
        System.out.println(" End test of query_disc2____________________________");
        
    }
    /**
     * tests GetNside method
     */
    public void testGetNside() {
        System.out.println(" Start test GetNside !!!!!!!!!!!!!!!!!!!!!");

        double pixsize = 0.3;
        long nside = pixTools.GetNSide(pixsize);
        System.out.println("Requared nside is "+nside);
        assertEquals("nside = " + nside, 1048576, nside, 1e-1);
        System.out.println(" End of GetNSide test _______________________");
    }
 
    /**
     *  test conversion of Ra Dec to polar coordinates
     */
    public void testRaDecToPolar() {
        System.out.println(" Start test RaDecToPolar !!!!!!!!!!!!!!!!!!!!!");
        double [] radec = new double[2];
        radec[0] = 312.115456;
        radec[1] = -1.153759;
        double[] polar = pixTools.RaDecToPolar(radec);
        assertEquals("theta = " + polar[0], 1.5909332201194137, polar[0], 1e-10);
        assertEquals("phi = " + polar[1], 5.447442353563491, polar[1], 1e-10);
        System.out.println("End test RaDecToPolar__________________________");
        
    }
}

