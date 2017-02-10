package gov.fnal.eag.healpix;

import javax.vecmath.Vector3d;

/**
 * Aggregates a vector and double[3][4] array which are the results of a
 * pixel center calculation.
 */
class PixInfo {
	final Vector3d pixVect;
	final double[][] pixVertex;

	PixInfo(Vector3d pixVect, double[][] pixVertex) {
		this.pixVect = pixVect;
		this.pixVertex = pixVertex;
	}
}
