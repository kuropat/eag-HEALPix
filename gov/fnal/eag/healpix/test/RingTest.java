package gov.fnal.eag.healpix.test;
import java.util.Collection;
import java.util.Iterator;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

import javax.vecmath.Vector3d;

import gov.fnal.eag.healpix.PixTools;

public class RingTest {
    private int nCheck;
    private int nFail;

 
    public Collection queryRing( int order, double ra, double dec,
                                 double radius ) {
        PixTools pixTools = PixTools.getInstance();
        long nside = orderToNside( order );
        double[] radec = new double[]{ra,dec};
        double[] polar = pixTools.RaDecToPolar(radec);
//        double theta = Math.PI / 2. - Math.toRadians(dec);
//        double phi = Math.toRadians(ra);
        Vector3d vec = pixTools.Ang2Vec( polar[0], polar[1] );
        return pixTools.query_ring( nside, vec, radius, 0, 0 );
    }

    public String report() {
        return "Checks: " + nCheck + "; "
             + "Failures: " + nFail;
    }

    static long orderToNside( int order ) {
        return 1L << order;
    }

    static Collection intersection( Collection c1, Collection c2 ) {
//    	System.out.println("n pixels 1="+c1.size()+" 2="+c2.size());
        Set s = new TreeSet( c1 );
        s.retainAll( c2 );
        return s;
    }

    /**
     * Haversine formula for spherical trigonometry.
     * This does not have the numerical instabilities of the cosine formula
     * at small angles.
     * <p>
     * This implementation derives from Bob Chamberlain's contribution
     * to the comp.infosystems.gis FAQ; he cites
     * R.W.Sinnott, "Virtues of the Haversine", Sky and Telescope vol.68,
     * no.2, 1984, p159.
     *
     * @param   ra1  right ascension of point 1 in radians
     * @param   dec1 declination of point 1 in radians
     * @param   ra2  right ascension of point 2 in radians
     * @param   dec2 declination of point 2 in radians
     * @return  angular separation of point 1 and point 2 in radians
     * @see  <http://www.census.gov/geo/www/gis-faq.txt>
     */
    static double haversineSeparationFormula( double ra1, double dec1,
                                              double ra2, double dec2 ) {
    	
        double sd2 = Math.sin( 0.5 * ( dec2 - dec1 ) );
        double sr2 = Math.sin( 0.5 * ( ra2 - ra1 ) );
        double a = sd2 * sd2 +
                   sr2 * sr2 * Math.cos( dec1 ) * Math.cos( dec2 );
        if ( Double.isNaN( a ) ) {
            return Double.NaN;
        }
        return a < 1.0 ? 2.0 * Math.asin( Math.sqrt( a ) )
                       : Math.PI;
    }
    /**
     * 
     * @param order
     * @param raDeg1 in degrees
     * @param decDeg1
     * @param raDeg2
     * @param decDeg2
     * @throws Exception
     */
    public void checkRing( int order, double raDeg1, double decDeg1,
            double raDeg2, double decDeg2 ) throws Exception {
    	double ra1 = raDeg1;
    	double dec1 = decDeg1 ;
    	double ra2 = raDeg2;
    	double dec2 = decDeg2;
    	double[] radec1 = {raDeg1,decDeg1};
    	double[] polar1 = new double[2];
    	double[] radec2 = {raDeg2,decDeg2};
    	double[] polar2 = new double[2];
    	long nside = orderToNside(order);
//    	double pixSize = Math.PI / (2.0 * nside); // in radians
    	long npix = 12*nside*nside;
    	double pixSize = Math.sqrt(4.0*Math.PI/npix);
    	System.out.println("pixSize="+pixSize+" nside="+nside);
//System.out.println("ra1="+ra1+" dec1="+dec1+" ra2="+ra2+" dec2="+dec2);
    	PixTools pixTools = PixTools.getInstance();
    	polar1 = pixTools.RaDecToPolar(radec1);
    	polar2 = pixTools.RaDecToPolar(radec2);
    
//System.out.println("thet1="+polar1[0]+" phi1="+polar1[1]+" thet2="+polar2[0]+" phi2="+polar2[1]);
    	Vector3d v1 = pixTools.Ang2Vec(polar1[0], polar1[1]); // center of 1
    	Vector3d v2 = pixTools.Ang2Vec(polar2[0], polar2[1]); // center of 2

// Distance between the two points.
    	double separation = haversineSeparationFormula( ra1, dec1, ra2, dec2 );
    	double separation2 = pixTools.AngDist(v1, v2);
    	System.out.println("separation1="+separation+" separation2="+separation2);
// Pick disc radius to be slightly larger than half the distance.
// The discs must therefore overlap, which means that the
// returned pixel lists must have at least one pixel in common.
    	double radius = separation * 0.5001;
    	double radius2 = separation2 * 0.5001;
    	radius = 20.0*pixSize;
    	System.out.println("radius="+radius+" radius2="+radius2);
    	Collection disc1 = queryRing( order, ra1, dec1, radius );
//    	Collection disc2 = queryRing( order, ra2, dec2, radius );
    	
    	Vector3d vP = new Vector3d(0.,0.,0.);
    	long cp = pixTools.vect2pix_ring(nside, v1);
    	v1 = pixTools.pix2vect_ring(nside, cp);
    	polar1 = pixTools.Vect2Ang(v1);
    	radec1 = pixTools.PolarToRaDec(polar1);
    	disc1 = queryRing( order, radec1[0], radec1[1], radius );
    	Iterator it = disc1.iterator(); 
    	long pix = 0;
    	System.out.println("number of pixels in ring "+disc1.size());
    	while (it.hasNext()) {
    		pix = ((Long)it.next()).longValue();
//    		polar1 = pixTools.pix2ang_ring(nside, pix);
//    		vP = pixTools.Ang2Vec(polar1[0], polar1[1]);
    		vP=pixTools.pix2vect_ring(nside, pix);
    		separation = pixTools.AngDist(v1, vP);
    		System.out.println("separation="+separation);
    	}
 
}

    public static void main( String[] args ) throws Exception {
        RingTest test = new RingTest();
        Random rnd = new Random( 2301L );
        double scale = 0.1; // degree
//        for ( int i = 0; i < 1; i++ ) {
//            double ra1 = rnd.nextDouble() * 360;
//            double dec1 = rnd.nextDouble() * 180 - 90;
//            double ra2 = Math.min( ra1 + rnd.nextDouble() * scale, 360 );
//            double dec2 = Math.max( Math.min( dec1 + rnd.nextDouble() * scale,
//                                              90 ), -90 );
//            for ( int order = 5; order < 12; order++ ) {
//                test.checkRing( order, ra1, dec1, ra2, dec2 );
//            }
//        }
//        test.checkRing( 8, 0., -89.983888, 180., -89.983888 );
    	test.checkRing( 10, 10.0, 0.0, 180., 0.1 );

        System.out.println( test.report() );
    }
    
}
