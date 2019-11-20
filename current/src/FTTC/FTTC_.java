package FTTC;

import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.*;
import ij.io.*;
import java.io.*;
import java.util.*;
import java.text.NumberFormat;
import Jama.Matrix;
import edu.emory.mathcs.jtransforms.fft.*;


/*
This is the main code for doing the FTTC calculation as an ImageJ plugin. 
The same package contains a plotting class as well. You will need it as well as the jama and jtransform to compile the binary. 
For more detail please see:
https://sites.google.com/site/qingzongtseng/tfm
or contact me. 
By TSENG Qingzong (qztseng /at/ gmail.com)
*/


public class FTTC_ implements PlugIn {

    private static double pixel;                             // pixel length in micron
    double dPixel;                                  //distance between each point of displacement (in pixel)
    private static double mu;                                   // the poission ratio of gel
    private static double E;                                       // Young's modulus of gel, in Pascal
    private static double lambda;                          //modified for test 11-02-16
    boolean drawPlot;
    private static int width, height;
    String dir, file;
    int magW, magH;

    public void run(String arg) {

        double[][] disXY = null;
        double[][] disX, disY, disXCF, disYCF, TractionXR, TractionYR;
        double[] Kx, Ky;                                     // the wavevector array
        double[] gridX, gridY;
        double D;
        int[] dim = null ;// the real distance between each point (in micron)
        Matrix TXY, G, H, G1, Gt, Ginv, dispXY, GtU;

        if (!getParams()) {
            return;
        }

        try {
            disXY = loadMatrixFromFile(dir + file, 4);
            //disXY is a 2D array with each row correspond to each vector, 4 columns correspond to x,y,dx,dy
        } catch (Exception e) {
            IJ.error(e.getMessage());
        }

        dim = plot_FTTC.getDimensions(disXY);
        //dim[0] will be the width(column)
        //dim[1] will be the height (row)
        //dim[3] indicates whether it is row major (1) or column major(0)
        magW = (dim[0]-1)*dim[2];
        magH = (dim[1]-1)*dim[2];
        if(width == 0 && height ==0){
            width = (dim[0]) * dim[2] + (int) (2 * disXY[0][0]);
            height = (dim[1]) * dim[2] + (int) (2 * disXY[0][1]);
        }

        disX = getSubMatrix(disXY, dim[1], dim[0], dim[3], 2);
        //disX are 2D arrays with disX[x][y] = "the displacement in X direction at vector(x,y)"
        disY = getSubMatrix(disXY, dim[1], dim[0], dim[3], 3);

        gridX = getCoord(disXY, dim[1], dim[0], dim[3], 0);
        gridY = getCoord(disXY, dim[1], dim[0], dim[3], 1);

        // convert distance in pixel into distance in meter (multiply by the pixel size)
        disX = scaleMatrix(disX, pixel);
        disY = scaleMatrix(disY, pixel);

        // get the distance between each data point (in pixel)
        dPixel = dim[2];

        //calculate the real distance between each data point (in meter)
        D = dPixel * pixel;

        // record the original matrix size before padding and adding imaginary part
        int nCol0 = disX[0].length;
        int nRow0 = disX.length;

        //start timer
        long startTime = System.currentTimeMillis();

        int nCol = nCol0;
        int nRow = nRow0;

        try {
            disX = paddingZero(disX, nCol * 2, nRow);         // padding 2 times in col for the imaginary part
            disY = paddingZero(disY, nCol * 2, nRow);
        } catch (Exception e) {
            IJ.error("padding failed");
        }

        //convert disX, disY into fourier space
        DoubleFFT_2D fftX = new DoubleFFT_2D(nRow, nCol);
        DoubleFFT_2D fftY = new DoubleFFT_2D(nRow, nCol);
        fftX.realForwardFull(disX);
        fftY.realForwardFull(disY);

        // suppress the net translation by setting the first 0th element in fourier space to zero
        disX[0][0] = 0;
        disX[0][1] = 0;
        disY[0][0] = 0;
        disY[0][1] = 0;

        //disXCF is a complex matrix now, the [0][0] is the real part of the first element, [0][1] is the imaginary part of the first element
        disXCF = disX;
        disYCF = disY;

        //setup the wavevector array
        Kx = new double[nCol];
        Ky = new double[nRow];

        for (int i = 0; i <= nCol / 2; i++) {
            Kx[i] = (2 * Math.PI / D) * ((double) i / (double) nCol);
        }
        for (int i = Math.round((float)nCol / 2 - 1); i > 0; i--) {              // the negative frequency part of the wavevector
            Kx[nCol - i] = -(2 * Math.PI / D) * ((double) i / (double) nCol);
        }
        for (int i = 0; i <= nRow / 2; i++) {
            Ky[i] = (2 * Math.PI / D) * ((double) i / (double) nRow);
        }
        for (int i = Math.round((float)nRow / 2 - 1); i > 0; i--) {              // the negative frequency part of the wavevector
            Ky[nRow - i] = -(2 * Math.PI / D) * ((double) i / (double) nRow);
        }

        //create the empty matrix to store the calculated traction force
        double[][] TractionXF = new double[nRow][nCol * 2];
        double[][] TractionYF = new double[nRow][nCol * 2];
        // Matrix to store temporarily the calculated traction force
        TXY = new Matrix(2, 1);
        // Matrix of the Green function
        G = new Matrix(2, 2);
        //Identity matrix for 0th order regularization
        H = Matrix.identity(2, 2);
        // for 0th order regularization this term will be constant. See the equation (11) from: Sabass Biophy. J 2008
        Double lambdaSq = lambda * lambda;
        H.timesEquals(lambdaSq);

        for (int j = 0; j < Ky.length; j++) {
            for (int i = 0; i < Kx.length; i++) {
                double k = Math.sqrt(Kx[i] * Kx[i] + Ky[j] * Ky[j]);
                if (i == (float)nCol / 2 + 1 || j == (float)nRow / 2 + 1) {
                    // at the nyquist frequency, set the off-diagonal element to zero, see Butler et al. Am J Physil Cell Physiol 2001
                    G.set(0, 1, 0);
                    G.set(1, 0, 0);
                } else if (i != 0 || j != 0) {
                    double gg = -mu * Kx[i] * Ky[j];
                    G.set(0, 1, -mu * Kx[i] * Ky[j]);
                    G.set(1, 0, -mu * Kx[i] * Ky[j]);
                }
                // the Green function
                double G0 = 2 * (1 + mu) / (E * Math.pow(k, 3));
                G.set(0, 0, (1 - mu) * (k * k) + mu * (Ky[j] * Ky[j]));
                G.set(1, 1, (1 - mu) * (k * k) + mu * (Kx[i] * Kx[i]));
                G.timesEquals(G0);

                Gt = G.transpose();
                G1 = Gt.times(G);
                G1 = G1.plus(H);
                Ginv = G1.inverse();

                //to get the displacement
                double[][] dd = new double[2][2];
                dd[0][0] = disXCF[j][i * 2];
                dd[0][1] = disXCF[j][i * 2 + 1];

                dd[1][0] = disYCF[j][i * 2];
                dd[1][1] = disYCF[j][i * 2 + 1];

                dispXY = new Matrix(dd);
                GtU = Gt.times(dispXY);
                TXY = Ginv.times(GtU);

                //Store the calculated traction
                TractionXF[j][i * 2] = TXY.get(0, 0);
                TractionXF[j][i * 2 + 1] = TXY.get(0, 1);
                TractionYF[j][i * 2] = TXY.get(1, 0);
                TractionYF[j][i * 2 + 1] = TXY.get(1, 1);
            }
        }
        //set the net traction(first element in the matrix) to zero
        TractionXF[0][0] = 0;
        TractionXF[0][1] = 0;
        TractionYF[0][0] = 0;
        TractionYF[0][1] = 0;

        //inverse FFT of traction matrix to get the force in real space
        fftX.complexInverse(TractionXF, true);
        fftY.complexInverse(TractionYF, true);

        // To get the real part of the calculated Force
        TractionXR = real(TractionXF);
        TractionYR = real(TractionYF);

        //stop timer
        long elapsedTime = System.currentTimeMillis() - startTime;

        IJ.showStatus("FTTC done in " + elapsedTime + "msec");
        IJ.log("FTTC done in " + elapsedTime + "msec");
        IJ.log("E (Young's modulus, in Pascal): " + E);
        IJ.log("lambda: " + lambda);
        IJ.log("pixel size(in micron): " + pixel * 1E6);

        /*draw the vector plot*/
        double[][] mag = calcMagnitude(TractionXR, TractionYR);
        Matrix mmag = new Matrix(mag);
        Matrix mmagt = mmag.transpose();
        double[][] magt = mmagt.getArray();
        if (drawPlot) {
            double[][] TracXY = saveFTTCdata(gridX, gridY, TractionXR, TractionYR, mag, dir + "Traction_" + file);
            //TracXY is 2D array organized like disXY, while TractionXR, Traction YR, mag are 2D arrays with array[x][y] correspond to vector(x,y)
            plot(TracXY, magt, file);

        } else {
            saveFTTCdata(gridX, gridY, TractionXR, TractionYR, mag, dir + "Traction_" + file);
        }

        saveFTTCparam(E, lambda, pixel, mu, dir + "FTTCparameters_" + file);


    }

    private double[][] calcMagnitude(double[][] TX, double[][] TY) {

        double[][] result = new double[TX.length][TX[0].length];

        for (int i = 0; i < TX.length; i++) {
            for (int j = 0; j < TX[0].length; j++) {
                result[i][j] = (float) Math.sqrt((TX[i][j] * TX[i][j]) + (TY[i][j] * TY[i][j]));
            }
        }

        return result;

    }

    private void plot(double[][] data, double mag[][], String title) {


        ImageProcessor ip = new ColorProcessor(width, height);
        int[] dim = plot_FTTC.getDimensions(data);
        double max = plot_FTTC.findMax2DArray(mag);
        double sc;
        plot_FTTC.colorMax = max;
        sc = 24/max;

        plot_FTTC.loadLut("S_Pet");
        plot_FTTC.drawVectors(ip, dim, data, mag, sc, plot_FTTC.colors);
        ImagePlus vp = new ImagePlus(title, ip);
        vp.show();
        plot_FTTC.makeScaleGraph(sc);

        //show the magnitude
        FloatProcessor ipm = new FloatProcessor(plot_FTTC.doubleArrayToFloat(mag));
        ipm.setInterpolationMethod(ImageProcessor.BICUBIC);
        ipm = (FloatProcessor) ipm.resize(magW, magH, false);
        ImageProcessor ipm2 = ipm.createProcessor(width, height);
        ipm2.setValue(0.0);
        ipm2.fill();
        ipm2.insert(ipm, (int)data[0][0], (int)data[0][1]);
        ImagePlus impM = new ImagePlus("Traction Magnitude_" + title, ipm2);
        impM.show();

        return;
    }

    private double[][] loadMatrixFromFile(String path, int col) throws Exception {

        ArrayList<String[]> aL = new ArrayList<String[]>(); 	// an array list to hold all vector info as String[]
        String[] cell;
        String line;
        int row = 0;
        double[][] matrix;
        NumberFormat nf = NumberFormat.getInstance();


        try {
            /* open the file */
            BufferedReader r = new BufferedReader(new FileReader(path));
            row = 0;
            while ((line = r.readLine()) != null) {
                line = line.trim();
                aL.add(line.split("\\s+"));     // each element of the array list: "line" is actually a string[4]
                row++;
            }
            r.close();
        } catch (Exception e) {
            throw new Exception("Unsupported file format.");
        }

        matrix = new double[row][col];
        Iterator<String[]> iter = aL.iterator();
        int counter = 0;

        /* go over all elements of aL*/
        while (iter.hasNext()) {
            cell = (String[]) iter.next();     // cell correspond to one element(string[4]) of aL, that is: one entry of data containing x,y,ux,uy, L
            //check if we have the fifth column which correspond to the vector magnitude in our data.
            if (cell.length<col){
                throw new Exception("Uncompatible data format");
            }else{
                for (int i = 0; i < col; i++) {
                    matrix[counter][i] = (Double)nf.parse(cell[i]).doubleValue();
                }
                counter++;
            }
        }

        return matrix;
    }

    static public double[][] paddingZero(double[][] matrix, int sizeX, int sizeY) throws Exception {

        if (matrix.length > sizeY || matrix[0].length > sizeX) {
            // IJ.log("exception!!");
            throw new Exception("matrix larger than specified size");

        }

        double[][] matP = new double[sizeY][sizeX];   //42x80
        for (int i = 0; i < matP.length; i++) {     //42
            for (int j = 0; j < matP[0].length; j++) {      //80
                if (i < matrix.length && j < matrix[0].length) {
                    matP[i][j] = matrix[i][j];
                } else {
                    matP[i][j] = 0.0;
                }

            }

        }

        return matP;
    }

    double[][] real(double[][] compMat) {

        double[][] realMat = new double[compMat.length][compMat[0].length / 2];

        for (int j = 0; j < realMat.length; j++) {
            for (int i = 0; i < realMat[0].length; i++) {
                realMat[j][i] = compMat[j][i * 2];
            }
        }

        return realMat;
    }

    private double[][] saveFTTCdata(double[] gridX, double[] gridY, double[][] tracX, double[][] tracY, double[][] mag, String path) {

        if (gridX.length != tracX[0].length || gridX.length != tracY[0].length || gridY.length != tracX.length || gridY.length != tracY.length) {
            IJ.error("gridX.length: " + gridX.length);
            IJ.error("gridY.length: " + gridY.length);
            IJ.error("tracX.length: " + tracX.length);
            IJ.error("tracX[0].length: " + tracX[0].length);
            IJ.error("tracY.length: " + tracY.length);
            IJ.error("tracY[0].length: " + tracY[0].length);
            IJ.error("the dimension of FTTC output is not the same as grid");
            return null;
        }
        double[][] TXY = new double[tracX.length * tracX[0].length][5];

        StringBuilder sb = new StringBuilder();

        for (int j = 0; j < tracX.length; j++) {
            for (int i = 0; i < tracX[0].length; i++) {
                TXY[tracX[0].length * j + i][0] = gridX[i];
                TXY[tracX[0].length * j + i][1] = gridY[j];
                TXY[tracX[0].length * j + i][2] = tracX[j][i];
                TXY[tracX[0].length * j + i][3] = tracY[j][i];
                TXY[tracX[0].length * j + i][4] = mag[j][i];

                sb.append(Double.toString(gridX[i]));
                sb.append(" ");
                sb.append(Double.toString(gridY[j]));
                sb.append(" ");
                sb.append(Double.toString(tracX[j][i]));
                sb.append(" ");
                sb.append(Double.toString(tracY[j][i]));
                sb.append(" ");
                sb.append(Double.toString(mag[j][i]));
                sb.append("\n");
            }
        }

        String s = sb.toString();
        IJ.saveString(s, path);
        return TXY;
    }

    private void saveFTTCparam(double _E, double _lambda, double _pixel, double _Poisson, String path) {

        StringBuilder sb = new StringBuilder();

        sb.append("Young's Modulus (E): ");
        sb.append(Double.toString(_E));
        sb.append("\n");
        sb.append("Regularization factor(lambda): ");
        sb.append(Double.toString(_lambda));
        sb.append("\n");
        sb.append("pixel size (in micron): ");
        sb.append(Double.toString(_pixel));
        sb.append("\n");
        sb.append("Poisson ratio: ");
        sb.append(Double.toString(_Poisson));
        sb.append("\n");


        String s = sb.toString();
        IJ.saveString(s, path);
        return;
    }

    private double[][] getSubMatrix(double[][] matrix, int row, int col, int rowMajor, int element) {

        double[][] subM = new double[row][col];
        for (int j = 0; j < row; j++) {
            for (int i = 0; i < col; i++) {
                if (rowMajor == 1) {
                    subM[j][i] = matrix[j * col + i][element];
                } else {
                    subM[j][i] = matrix[i * row + j][element];
                }
            }
        }

        return subM;
    }

    private double[] getCoord(double[][] matrix, int row, int col, int rowMajor, int element) {

        double[] Coord;

        if (element == 0) {    // getting x coordinates
            Coord = new double[col];
            for (int i = 0; i < col; i++) {
                if (rowMajor == 1) {   //y was fixed first
                    Coord[i] = matrix[i][0];
                } else {
                    Coord[i] = matrix[i * row][0];
                }
            }
        } else if (element == 1) {      //getting y coordinates
            Coord = new double[row];
            for (int i = 0; i < row; i++) {
                if (rowMajor == 1) {   //y was fixed first
                    Coord[i] = matrix[i * col][1];
                } else {
                    Coord[i] = matrix[i][1];
                }
            }
        } else {
            IJ.error("error!");
            Coord = new double[0];
        }

        return Coord;

    }

    private double[][] scaleMatrix(double[][] mat, double scale) {

        double[][] sMat = new double[mat.length][mat[0].length];
        for (int j = 0; j < mat.length; j++) {
            for (int i = 0; i < mat[0].length; i++) {
                sMat[j][i] = mat[j][i] * scale;
            }
        }

        return sMat;
    }

    private boolean getParams() {

        GenericDialog gd = new GenericDialog("FTTC");
        if (pixel == 0.0D) {
            pixel = 0.09E-6;
        }
        gd.addNumericField("pixel size(in micron)", pixel * 1E6, 3);
        if (mu == 0.0D) {
            mu = 0.5;
        }
        gd.addNumericField("poisson ratio of the gel", mu, 1);
        if (E == 0.0D) {
            E = 5000;
        }
        gd.addNumericField("Young's modulus of the gel(in Pascal)", E, 2);

        gd.addNumericField("regularization factor", lambda, 15, 15, "");
        gd.addCheckbox("plot result?", true);
        if (width == 0 || height == 0) {
            width = 1000;
            height = 1000;
        }
        gd.addNumericField("plot width:", width, 0);
        gd.addNumericField("plot height:", height, 0);

        gd.showDialog();

        pixel = (1E-6) * gd.getNextNumber();
        mu = gd.getNextNumber();
        E = gd.getNextNumber();
        lambda = gd.getNextNumber();
        drawPlot = gd.getNextBoolean();
        width = (int) gd.getNextNumber();
        height = (int) gd.getNextNumber();

        if (gd.wasCanceled()) {
            return false;
        }

        OpenDialog od = new OpenDialog("Select the displacement file", "");
        if (od.getDirectory() == null || od.getFileName() == null) {
            return false;
        }
        dir = od.getDirectory();
        file = od.getFileName();

        return true;
    }

       static public StringBuffer generatePIVToPrint(double[][] data) {
        NumberFormat nf = NumberFormat.getInstance();
        nf.setMaximumFractionDigits(12);
        nf.setMinimumFractionDigits(12);

        StringBuffer info = new StringBuffer();

        for (int i = 0; i < data.length; i++) {
            for (int c = 0; c < data[0].length; c++) {
                info.append(nf.format(data[i][c]));
                info.append(" ");
            }
            info.append("\n");
        }

        return info;
    }


    static public StringBuffer generateArrayToPrint(double[] data) {
        NumberFormat nf = NumberFormat.getInstance();
        nf.setMaximumFractionDigits(12);
        nf.setMinimumFractionDigits(12);

        StringBuffer info = new StringBuffer();

        for (int i = 0; i < data.length; i++) {
            info.append(nf.format(data[i]));
            info.append("\n");
        }

        return info;
    }
}
