package gov.fnal.eag.healpix.test;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;

import gov.fnal.eag.healpix.BitManipulation;
import gov.fnal.eag.healpix.PixTools;
import javax.vecmath.Vector3d;
import junit.framework.*;

public class QueryDiscTest extends TestCase {
	public void testQueryDisc () {
	long nside = 32;
	
	int nest = 0;
	int inclusive = 0;
	double radius = Math.PI;
	double radius1 = Math.PI/2.;
    PixTools pt = new PixTools();
    int npix = (int) pt.Nside2Npix(nside);
    double res = pt.PixRes(nside); // pixel size in radians
    System.out.println("res="+res);
    double pixSize = Math.toRadians(res/3600.0); // pixel size in radians
    System.out.println("pixSize="+pixSize+" rad");

    
    ArrayList fullSky = pt.query_disc(nside, new Vector3d(0., 0., 1.), radius, nest, inclusive);
    ArrayList firstHalfSky = pt.query_disc(nside, new Vector3d(0., 0., 1.), radius1, nest, inclusive);
    ArrayList secondHalfSky = pt.query_disc(nside, new Vector3d(0., 0., -1.), radius1, nest, inclusive);
    firstHalfSky.addAll(secondHalfSky);
    HashSet pixHalfsUnique = new HashSet(firstHalfSky);
    ArrayList pixHalfsList = new ArrayList(pixHalfsUnique);
    Collections.sort(pixHalfsList);
    Collections.sort(fullSky);

    long listL = Math.min(fullSky.size(),pixHalfsList.size() );
    assertEquals(npix,fullSky.size());
    assertEquals(npix,listL);
    for ( int i=0; i< listL; i++) {

    assertEquals(fullSky.get(i),pixHalfsList.get(i));
    }
    
    double[] ang_tup = { 0., 0. };
    Vector3d v = new Vector3d(1., 0., 0.);
    Vector3d v1 = new Vector3d(-1., 0., 0.);
	ang_tup = pt.Vect2Ang(v);

   firstHalfSky = pt.query_disc(nside, new Vector3d(1., 0., 0.), radius1, nest, inclusive);
   secondHalfSky = pt.query_disc(nside, new Vector3d(-1., 0., 0.),radius1, nest, inclusive);
    firstHalfSky.addAll(secondHalfSky);
    pixHalfsUnique = new HashSet(firstHalfSky);
    pixHalfsList = new ArrayList(pixHalfsUnique);
    
    Collections.sort(pixHalfsList);
//    System.out.println("full size="+fullSky.size()+" half size="+pixHalfsList.size());
    listL = Math.min(fullSky.size(),pixHalfsList.size() );
    assertEquals(npix,fullSky.size());
    assertEquals(npix,listL);
    for ( int i=0; i< listL; i++) {
//        System.out.println( "i="+i+" "+fullSky.get(i)+" "+pixHalfsList.get(i));
        assertEquals(fullSky.get(i),pixHalfsList.get(i));
        }


    firstHalfSky = pt.query_disc(nside, new Vector3d(0., 1., 0.), radius1, nest, inclusive);
    secondHalfSky = pt.query_disc(nside, new Vector3d(0., -1., 0.), radius1, nest, inclusive);
    firstHalfSky.addAll(secondHalfSky);
    pixHalfsUnique = new HashSet(firstHalfSky);
    pixHalfsList = new ArrayList(pixHalfsUnique);
    Collections.sort(pixHalfsList);
//    System.out.println("full size="+fullSky.size()+" half size="+pixHalfsList.size());
    listL = Math.min(fullSky.size(),pixHalfsList.size() );
    assertEquals(npix,fullSky.size());

    for ( int i=0; i< listL; i++) {

        assertEquals(fullSky.get(i),pixHalfsList.get(i));
        }
        
}
	
	   public void testQueryCircle(){
	       PixTools pt = new PixTools();

	       long nside = 512;

	       double angle = Math.toRadians(0.011451621372724687);
//	       double angle = Math.toRadians(0.11451621372724687); 
	       Vector3d v = new Vector3d(0.8956388362603873,
	   -1.838600645782914E-4, 0.44478201534866);

	       //convert vector to IPIX
	       long ipix = pt.vect2pix_ring(nside, v);
	       //and query circle
	       ArrayList r = pt.query_disc(nside,v,angle,0,1);

	       //now make test that IPIX is in Circle, this will fail
	       System.out.println("ipix = "+ipix);
	       System.out.println("disc: "+r);
	   //ipix = 875520
	   //disc: [873471, 871424, 875519, 873472]

	       assertTrue("pixel not found in disc",r.contains(new Long(ipix)));
	       ArrayList r1 = pt.query_disc(nside,v,angle,0,0);
	       assertTrue("pixel not found in disc",r1.contains(new Long(ipix)));
	      }
	    /**
	     * additional tests of Query_Disc method
	     * @throws Exception
	     */
	    public void testQuery_Disc_1() throws Exception {
	        System.out.println(" Start test query_disc_1 !!!!!!!!!!!!!!!!!!!!!");
	        long[] case1 = {282, 271, 280, 269};
	        long[] case2 = { 295, 293, 292};
	        long[] case3 = {312, 301, 306, 305, 300, 295, 304, 294, 293, 282, 291, 292, 271, 289, 270, 269, 267};
	        long[] case4 = { 63, 127, 191, 255, 62, 61, 126, 125, 190, 189, 254, 253};

	        PixTools pt = new PixTools();
	        long nside = 8;
	        int nest = 1;
	        long ipix = 0;
	        System.out.println("case 1");
	        double pixres= Math.toRadians(pt.PixRes(nside)/3600.);
	        System.out.println("pix res="+pixres+" rad");

	        int inclusive = 1;
	        double radius = Math.toRadians(1.932);

	        double[] radec = new double[]{5.64001,-4.57566};
	        double[] polar = pt.RaDecToPolar(radec);
	        Vector3d v = pt.Ang2Vec(polar[0], polar[1]);
	        ArrayList pixlist;
	        pixlist = pt.query_disc(nside, v, radius, nest, inclusive);

	        int nlist = pixlist.size();
	        long pixel =pt.vect2pix_nest(nside,v);
	        
	        System.out.println(" radius="+radius+" center pixel="+pixel);
	        for (int i = 0; i < nlist; i++) {
	            ipix = ((Long) pixlist.get(i)).longValue();
	          System.out.println("i="+i+" pixel="+ipix);
	            assertEquals("pixel = " + ipix, case1[i], ipix, 1e-10);
	        }
	        System.out.println("case 2");
	           radec = new double[]{351.04536,1.74210};
	           polar = pt.RaDecToPolar(radec);
	           v = pt.Ang2Vec(polar[0], polar[1]);
	           radius = Math.toRadians(1.932);

	           pixlist = pt.query_disc(nside, v, radius, nest, inclusive);

	           nlist = pixlist.size();
	           pixel =pt.vect2pix_nest(nside,v);
	            System.out.println(" radius="+radius+" center pixel="+pixel);
	            for (int i = 0; i < nlist; i++) {
	                ipix = ((Long) pixlist.get(i)).longValue();
	              System.out.println("i="+i+" pixel="+ipix);
	              assertEquals("pixel = " + ipix, case2[i], ipix, 1e-10);
	            }
	            System.out.println("case 3");
	            radius = Math.toRadians(10.68);
	            radec = new double[]{348.86174 ,-0.07390};
	               polar = pt.RaDecToPolar(radec);
	               v = pt.Ang2Vec(polar[0], polar[1]);
	               pixlist = pt.query_disc(nside, v, radius, nest, inclusive);

	               nlist = pixlist.size();
	               pixel =pt.vect2pix_nest(nside,v);
	                System.out.println(" radius="+radius+" center pixel="+pixel);
	                for (int i = 0; i < nlist; i++) {
	                    ipix = ((Long) pixlist.get(i)).longValue();
	                  System.out.println("i="+i+" pixel="+ipix);
	                assertEquals("pixel = " + ipix, case3[i], ipix, 1e-10);
	                } 
	                //
	                System.out.println(" case 4");
	                radius = Math.toRadians(4.134);

	                radec = new double[]{182.95228, 89.43585};
	                   polar = pt.RaDecToPolar(radec);
	                   v = pt.Ang2Vec(polar[0], polar[1]);
	                   pixlist = pt.query_disc(nside, v, radius, nest, inclusive);

	                   nlist = pixlist.size();
	                   pixel =pt.vect2pix_nest(nside,v);
	                    System.out.println(" radius="+radius+" center pixel="+pixel);
	                    for (int i = 0; i < nlist; i++) {
	                        ipix = ((Long) pixlist.get(i)).longValue();
	                      System.out.println("i="+i+" pixel="+ipix);
	                  assertEquals("pixel = " + ipix, case4[i], ipix, 1e-10);
	                    }  
	        

	        System.out.println(" End query_disk_1 test ______________________________");
	    }
}
