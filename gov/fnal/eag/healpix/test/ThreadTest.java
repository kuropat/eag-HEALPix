package gov.fnal.eag.healpix.test;


import gov.fnal.eag.healpix.PixTools;
import java.util.List;
import java.util.Random;
import javax.vecmath.Vector3d;
/**
 * The thread proof test by Mark Teylor
 * @author Mark Taylor
 *
 */
public class ThreadTest {
    private static PixTools pixTools = new PixTools();
    public void testPixtools() throws InterruptedException {
        workPixtools( false );
        workPixtools( true );
    }

    private void workPixtools( boolean useThreads )
            throws InterruptedException {
        int nw = 40;
        Worker[] workers = new Worker[ nw ];
        Random rnd = new Random( 23001L );
        int nside = 6 + rnd.nextInt( 3 );
        Vector3d vector = new Vector3d( rnd.nextDouble(),
                                        rnd.nextDouble(),
                                        rnd.nextDouble() );
        double radius = 0.1 * rnd.nextDouble();
        int nest = 1;
        int inclusive = 1;
        for ( int iw = 0; iw < nw; iw++ ) {
            workers[ iw ] =
                new Worker( nside, vector, radius, nest, inclusive );
        }
        if ( useThreads ) {
            for ( int iw = 0; iw < nw; iw++ ) {
                workers[ iw ].start();
            }
            for ( int iw = 0; iw < nw; iw++ ) {
                workers[ iw ].join();
            }
        }
        else {
            for ( int iw = 0; iw < nw; iw++ ) {
                workers[ iw ].run();
            }
        }
        List result0 = workers[ 0 ].result_;
        for ( int iw = 1; iw < nw; iw++ ) {
            assertEquals( "Fail at " + iw + " for useThreads=" + useThreads,
                          result0, workers[ iw ].result_ );
        }
    }

    private static void assertEquals( String msg, List l1, List l2 ) {
        if ( ! l1.equals( l2 ) ) {
            throw new AssertionError( msg + "(" + l1 + " != " + l2 + ")" );
        }
    }

    private static class Worker extends Thread {
        private final long nside_;
        private final Vector3d vector_;
        private final double radius_;
        private final int nest_;
        private final int inclusive_;
        private List result_;

        Worker( long nside, Vector3d vector, double radius,
                int nest, int inclusive ) {
            nside_ = nside;
            vector_ = vector;
            radius_ = radius;
            nest_ = nest;
            inclusive_ = inclusive;
        }

        public void run() {

            result_ = pixTools.query_disc( nside_, vector_, radius_, nest_, inclusive_ );
        }
    }

    public static void main( String[] args ) throws InterruptedException {
        new ThreadTest().testPixtools();
        System.out.println(" if no interrupts the test is complete successfully");
    }
}


