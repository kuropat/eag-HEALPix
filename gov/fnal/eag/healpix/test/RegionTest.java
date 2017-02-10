package gov.fnal.eag.healpix.test;

import java.util.ArrayList;
import javax.vecmath.Vector3d;
import gov.fnal.eag.healpix.PixTools;
import gov.fnal.eag.healpix.Region;
import gov.fnal.eag.healpix.BitManipulation;
import junit.framework.TestCase;

/**
 * @author kuropat
 * Junit methods to test Region class
 */
public class RegionTest extends TestCase {
	private final static double TWOPI = 2.0 * Math.PI;
	private final static double epsilon = 1.0e-10;
	private BitManipulation bm = new BitManipulation();
	/**
	 * test default constructor
	 */
	public void testRegion() {
		double xMin = 0.;
		double xMax = 30.;
		double yMin = 0.;
		double yMax = 10.;
		PixTools pt = new PixTools();
		Region rg = new Region(xMin,xMax,yMin,yMax);
		assertNotNull(rg);
		double x = -10.;
		double y = 15.;
		assertFalse(rg.inReg(x,y));
		x = 10.;
		y = 5.;
		assertTrue(rg.inReg(x,y));
		xMax = -10.;
		ArrayList vert = rg.getVertices();
		double[][] vertPol = rg.getPolReg();
		for ( int ind=0; ind<vert.size(); ind++) {
			Vector3d vv = (Vector3d) vert.get(ind);
			double [] vvAng = pt.Vect2Ang(vv);
			if (vertPol[ind][1] < 0) vertPol[ind][1] += TWOPI;
			double comp =  bm.MODULO(vvAng[1], TWOPI) - epsilon;
			assertEquals("theta ",vertPol[ind][0],vvAng[0], 1.0e-5);
			assertEquals("phi ="+vertPol[ind][1],vertPol[ind][1],vvAng[1], 1.0e-5);
		}
		xMin = 20.;
		xMax = 95.;
		Region rg1 = new Region(xMin,xMax,yMin,yMax);
		assertNotNull(rg1);
		x = 45.;
		y = 5.;
		assertTrue(rg1.inReg(x,y));
	}
	/**
	 * test pixelization
	 * 
	 */
	public void testPixelize() {
		System.out.println("test pixelize");
		ArrayList pixels = new ArrayList();
		double xMin = 10.;
		double xMax = 60.;
		double yMin = -20.0;
		double yMax = 0.;
		PixTools pt = new PixTools();
		Region rg = new Region(xMin,xMax,yMin,yMax);
		double[][] regCoord = rg.getPolReg();
		for (int i = 0; i<regCoord.length; i++ ) {
			System.out.println("thet="+regCoord[i][0]+" phi="+regCoord[i][1]);
		}
		double resolution = 10.*60.*60.; // resolution in arcsec (= 10 degrees)
		
		try {
			pixels = rg.pixelize(resolution);
			long nside = pt.GetNSide(resolution);
			int npix = pixels.size();
			assertFalse(npix == 0);
			System.out.println("npix="+npix);
			for (int i=0; i<npix; i++) {
				long pix = ((Long) pixels.get(i)).longValue();
				System.out.println("pixel="+pix);
				double[] pixpos = pt.pix2ang_ring(nside,pix);
				System.out.println("theta="+pixpos[0]+" phi="+pixpos[1]);
				double[] radec = pt.PolarToRaDec(pixpos);
				double[][] pixvert = pt.pix2vertex_ring(nside,pix);
				System.out.println("corners");
				for (int j=0; j<pixvert[0].length; j++) {
					double x = pixvert[0][j];
					double y = pixvert[1][j];
					double z = pixvert[2][j];
					double[] pol = pt.xyzToPolar(x,y,z);
					double[] radec1 = pt.PolarToRaDec(pol);
					System.out.println("ra= "+radec1[0]+" dec="+radec1[1]);
				}
				System.out.println();
				
			}
		
		} catch (Exception e) {
			System.err.println("Exception in pixelize");
			e.printStackTrace();
		}

	}
}
