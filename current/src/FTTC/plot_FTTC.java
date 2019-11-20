package FTTC;

import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.*;
import ij.io.OpenDialog;
import java.util.*;
import java.io.*;
import java.awt.image.*;
import java.text.*;


/*
This is the code for plotting the FTTC result. 
For more detail please see:
https://sites.google.com/site/qingzongtseng/tfm
or contact me. 
By TSENG Qingzong (qztseng /at/ gmail.com)
*/

public class plot_FTTC implements PlugIn {

    double[][] data;
    double[][] magArr;
    double scale = 0D;
    static double colorMax = 600D;
    static double max;
    ImagePlus vp, vpm;
    ImageProcessor ip, ipm, ipm2;
    String lut = "S_Pet";
    static String[] colors;
    static String[] spet = {"4 1 4","4 1 4","7 1 8","9 1 13","10 1 17","13 1 21","16 1 25","18 1 30","19 1 34","22 1 38","25 1 42","27 1 47","28 1 51","31 1 55","34 1 59","35 1 64","38 1 67","41 1 71","43 1 75","44 1 80","47 1 84","50 1 88","52 1 92","53 1 97","56 1 101","59 1 105","61 1 109","62 1 114","64 1 118","67 1 122","68 1 127","68 1 127","67 1 131","64 1 135","62 1 139","61 1 143","59 1 147","56 1 151","53 1 155","52 1 159","50 1 163","47 1 167","44 1 171","43 1 175","41 1 179","38 1 183","35 1 187","34 1 191","31 1 194","28 1 198","27 1 202","25 1 206","22 1 210","19 1 214","18 1 218","16 1 222","13 1 226","10 1 230","9 1 234","7 1 238","4 1 242","2 1 246","1 1 254","1 1 254","1 5 246","1 10 242","1 14 238","1 19 234","1 23 230","1 28 226","1 32 222","1 37 218","1 41 214","1 46 210","1 50 206","1 55 202","1 59 198","1 64 194","1 67 191","1 72 187","1 76 183","1 81 179","1 85 175","1 90 171","1 94 167","1 99 163","1 103 159","1 108 155","1 112 151","1 117 147","1 121 143","1 126 139","1 130 135","1 135 131","1 139 127","1 144 122","1 148 118","1 153 114","1 157 109","1 162 105","1 166 101","1 171 97","1 175 92","1 180 88","1 184 84","1 189 80","1 192 75","1 203 71","1 203 67","10 205 64","19 207 59","28 209 55","37 211 51","46 213 47","55 215 42","64 217 38","72 219 34","81 221 30","90 223 25","99 225 21","108 227 17","117 229 13","126 231 11","135 233 8","144 235 5","153 237 2","162 239 1","171 241 1","180 242 1","189 244 1","197 245 1","206 247 1","215 248 1","224 250 1","233 251 1","242 253 1","254 254 1","254 254 1","254 253 1","254 251 1","254 249 1","254 247 1","254 245 1","254 243 1","254 241 1","254 239 1","254 237 1","254 235 1","254 233 1","254 231 1","254 229 1","254 227 1","254 225 1","254 223 1","254 221 1","254 219 1","254 217 1","254 215 1","254 213 1","254 211 1","254 209 1","254 207 1","254 205 1","254 203 1","254 201 1","254 199 1","254 197 1","254 195 1","254 193 1","254 191 1","254 190 1","254 188 1","254 186 1","254 184 1","254 182 1","254 180 1","254 178 1","254 176 1","254 174 1","254 172 1","254 170 1","254 168 1","254 166 1","254 164 1","254 162 1","254 160 1","254 158 1","254 156 1","254 154 1","254 152 1","254 150 1","254 148 1","254 146 1","254 144 1","254 142 1","254 140 1","254 138 1","254 136 1","254 134 1","254 132 1","254 130 1","254 128 1","254 126 1","254 124 1","254 122 1","254 120 1","254 118 1","254 116 1","254 114 1","254 112 1","254 110 1","254 108 1","254 106 1","254 104 1","254 102 1","254 100 1","254 98 1","254 96 1","254 94 1","254 92 1","254 90 1","254 88 1","254 86 1","254 84 1","254 82 1","254 80 1","254 78 1","254 76 1","254 74 1","254 72 1","254 69 1","254 66 1","254 64 1","254 62 1","254 59 1","254 56 1","254 54 1","254 51 1","254 48 1","254 45 1","254 43 1","254 40 1","254 37 1","254 34 1","254 32 1","254 29 1","254 26 1","254 23 1","254 21 1","254 18 1","254 15 1","254 12 1","254 10 1","254 7 1","254 4 1","254 4 1"};
    static byte[] lutR, lutG, lutB;
    static IndexColorModel cm;
    boolean mag = false;
    boolean autoScale = false;
    int plotW = 0, plotH = 0;
    int pivW, pivH;
    boolean scalebar = false;
    int[] dimensions;		//record data dimensions as array (nColumn, nRow, vectorSpace)
    int maxLenPixel = 24;       //The longest arrow size in pixel in autoscaled plot

    
    public void run(String arg) {

        String[] pathname;

        try {
            pathname = getFilePath();
        } catch (Exception e) {
            return;
        }

        if (!getParams()) {
            return;
        }

        try {
            data = load2DArrayFromFile(pathname[0] + pathname[1]);
        } catch (Exception e) {
            IJ.error(e.getMessage());
        }

        /*load the lut*/
        loadLut(lut);

        /*get the data dimentsion and set the plot dimension*/
        dimensions = getDimensions(data);
        //dimensions[0] is the number of data points in each row (nx)
        //dimensions[1] is the number of data points in each column (ny)
        //dimensions[2] is the spacing between each data point (in pixel)

        /*Calculate the piv width and height, as well as the plot width and height*/
        pivW = (dimensions[0]-1)*dimensions[2];
        pivH = (dimensions[1]-1)*dimensions[2];
        if (plotH == 0 && plotW == 0) {
            plotW = (dimensions[0]) * dimensions[2] + (int) (2 * data[0][0]);
            plotH = (dimensions[1]) * dimensions[2] + (int) (2 * data[0][1]);
        }

        /*calculate the magnitude(norm of the vector)*/
        magArr = get2DElement(data, dimensions[0], dimensions[1], 4);

        if (autoScale || scalebar) {
            max = findMax2DArray(magArr);
            if (autoScale) {
                scale = maxLenPixel / max;
                colorMax = max;
            }
        }

        /*draw the magnitude plot*/
        if (mag) {
            ipm = new FloatProcessor(doubleArrayToFloat(magArr));
            ipm.setInterpolationMethod(ImageProcessor.BICUBIC);
            ipm = (FloatProcessor) ipm.resize(pivW, pivH, false);
            ipm.setColorModel(cm);
            ipm2 = ipm.createProcessor(plotW, plotH);
            ipm2.setValue(0.0);
            ipm2.fill();
            ipm2.insert(ipm, (int)data[0][0], (int)data[0][1]);
            vpm = new ImagePlus("Magnitude map_" + pathname[1], ipm2);
            vpm.show();
        }

        /*draw the vector plot*/
        ip = new ColorProcessor(plotW, plotH);
        drawVectors(ip, dimensions, data, magArr, scale, colors);
        vp = new ImagePlus("Vector plot_" + pathname[1], ip);
        vp.show();

        /*draw the scalebar*/
        if (scalebar) {
            makeScaleGraph(scale);
        }

        return;
    }

    protected static void makeScaleGraph(double sc) {

        ImagePlus impS;
        ColorProcessor cpScale;
        int nScaleArrows = 5;
        int colorBarHeight = 500;
        int colorBarWidth = 30;
        int sOffX = 10;
        int sOffY = 20;
        int scaleGraphWidth;
        scaleGraphWidth = (int) (colorMax * sc) + sOffX * 2;

        cpScale = makeVerticalColorBar(colorBarWidth, colorBarHeight, scaleGraphWidth, sOffY);
        addScaleArrows(cpScale, nScaleArrows, colorBarHeight, sc, sOffX, sOffY);
        impS = new ImagePlus("Scale Graph", cpScale);
        impS.show();

        return;
    }

    private static void addScaleArrows(ColorProcessor cp, int nArrows, int bHeight, double vs, int sOffX, int sOffY) {


        double[][] sArr = new double[nArrows + 1][4];
        double[][] sMagArr = new double[1][nArrows + 1];
        int[] dim = {1, nArrows, bHeight / nArrows};

        for (int i = 0; i <= nArrows; i++) {
            sArr[i][0] = sOffX;
            sArr[i][1] = sOffY + (bHeight / nArrows) * i;
            sArr[i][3] = 0;
            sArr[i][2] = colorMax * ((nArrows - i) / (double) nArrows);
            sMagArr[0][i] = sArr[i][2];

        }
        drawVectors(cp, dim, sArr, sMagArr, vs, colors);

        Color c = Color.white;
        cp.setColor(c);
        Font font = new Font("SansSerif", Font.BOLD, 12);
        cp.setFont(font);
        cp.setAntialiasedText(true);

        for (int i = 0; i <= nArrows; i++) {
            int xLabel = sOffX;
            int yLabel = (int) sArr[i][1] - 3;
            String label;
            label = IJ.d2s(colorMax * ((nArrows - i) / (double) nArrows), 0);

            cp.drawString(label, xLabel, yLabel);

        }

    }

    private static ColorProcessor makeVerticalColorBar(int _width, int _height, int offX, int offY) {

        //create a 0-255 vertical ramp image. Modified from the ImageJ newImage.class
        byte[] pixels = new byte[_width * _height];

        byte[] ramp = new byte[_height];
        int k = 0;
        for (int i = _height - 1; i >= 0; i--) {
            ramp[k] = (byte) ((i * 256.0) / _height);
            k++;

        }
        int offset;
        int a = 0;

        for (int y = 0; y < _height; y++) {
            offset = y * _width;
            for (int x = 0; x < _width; x++) {
                pixels[offset++] = ramp[a];
            }
            a++;
        }

        ImageProcessor bp = new ByteProcessor(_width, _height, pixels, null);
        bp.setColorModel(cm);
        ImageProcessor bp2 = bp.createProcessor(_width + offX, _height + offY);
        bp2.setValue(0.0);
        bp2.fill();
        bp2.insert(bp, offX, offY);
        ImagePlus temp = new ImagePlus("test", bp2);
        ColorProcessor bp2c = (ColorProcessor) temp.getProcessor().convertToRGB();
        temp.changes = false;
        temp.close();

        return bp2c;
    }

    private String[] getFilePath() throws Exception {
        String[] pn = new String[2];

        OpenDialog od = new OpenDialog("Select the data for plotting", "");
        if (od.getDirectory() == null || od.getFileName() == null) {
            throw new Exception("No file selected");
        }
        pn[0] = od.getDirectory();
        pn[1] = od.getFileName();

        return pn;
    }
//
//    double[][] readArrayFromFile(String pathname)
//            throws IOException {
//
//        int numOfRows = 0;
//        int i = 0;
//        int numOfCol = 0;
//        int j = 0;
//        double[] d = new double[25];
//        // array list for conveniently appending the parsed rows
//        // do not use a LinkedList here: bad performance
//        ArrayList al = new ArrayList();
//        // count lines that contain Tecplot key words
//        BufferedReader br = new BufferedReader(new FileReader(pathname));
//
//        // configuring StreamTokenizer
//        StreamTokenizer st = new StreamTokenizer(br);
//        st.eolIsSignificant(true);
//        st.resetSyntax();
//        st.wordChars('0', '9');
//        st.wordChars('-', '-');
//        st.wordChars('+', '+');
//        st.wordChars('e', 'e');
//        st.wordChars('E', 'E');
//        st.wordChars('.', '.');
//        st.whitespaceChars(' ', ' ');
//        st.whitespaceChars('\t', '\t');
//        int type = -1;
//        while ((type = st.nextToken()) != StreamTokenizer.TT_EOF) {
//            switch (type) {
//                case StreamTokenizer.TT_WORD:
//                    d[j] = Double.parseDouble(st.sval);
//                    j++;
//                    break;
//                case StreamTokenizer.TT_EOL:
//                    // at the end of the line, the data is appended at the ArrayList
//                    // use clone() to copy the object and not just its reference
//                    al.add(data.clone());
//                    numOfRows++;
//                    numOfCol = j;
//                    j = 0;
//                    break;
//                default:
//                    break;
//            }
//        }
//        br.close();
//        // copying the data from the ArrayList into a double-array
//        double[][] array = new double[numOfRows][numOfCol];
//        for (i = 0; i < numOfRows; ++i) {
//            d = (double[]) al.get(i);
//            System.arraycopy(data, 0, array[i], 0, numOfCol);
//        }
//        al.clear();
//        return (array);
//    }

    private boolean getParams() {

        GenericDialog gd = new GenericDialog("Vector plot");
        gd.addMessage("vector plot parameters:");
        gd.addCheckbox("Autoscale vector plot?", true);
        gd.addMessage("Set the following two scale factors if autoscale was not set:");
        gd.addNumericField("vector_scale", 1, 0);
        gd.addNumericField("max vector value", 500, 1);
        gd.addMessage("------------------------");
        gd.addNumericField("plot_width (pixel):", 0, 0);
        gd.addNumericField("plot_height (pixel):", 0, 0);
        gd.addMessage("(Plot dimension will be approximated the from the data\nif both width and height are set to 0)");
        gd.addMessage("------------------------");
        gd.addCheckbox("show scale graph?", true);
        gd.addCheckbox("draw magnitude(vector norm) plot?", false);
        gd.addStringField("LUT for color coding: ", "S_Pet");
        gd.showDialog();

        autoScale = gd.getNextBoolean();
        scale = (double) gd.getNextNumber();
        colorMax = (double) gd.getNextNumber();
        plotW = (int) gd.getNextNumber();
        plotH = (int) gd.getNextNumber();
        scalebar = gd.getNextBoolean();
        mag = gd.getNextBoolean();
        lut = gd.getNextString();
        
        if (gd.wasCanceled()) {
            return false;
        }

        return true;
    }

    protected static int[] getDimensions(double[][] array) {
        int[] dm = new int[4];                               //dm[2] is the spacing of vector (points), dm[0] is the number of points in x, dm[1] is np in y

        if ((array[1][0] - array[0][0]) == 0) {				// x was fixed first
            dm[3] = 0; // column major
            dm[2] = (int) (array[1][1] - array[0][1]);
            for (int i = 1; i < array.length; i++) {
                if (array[i][0] == array[i - 1][0]) {
                    dm[1]++;
                } else {
                    break;
                }
            }
            dm[1]++;
            dm[0] = array.length / dm[1];
        } else {
            dm[3] = 1; // row major
            dm[2] = (int) (array[1][0] - array[0][0]);
            for (int i = 1; i < array.length; i++) {
                if (array[i][1] == array[i - 1][1]) {
                    dm[0]++;
                } else {
                    break;
                }
            }
            dm[0]++;
            dm[1] = array.length / dm[0];
        }
        return dm;
    }

    protected static void drawVectors(ImageProcessor ip, int[] dim, double[][] d, double[][] magArr, double vs, String[] lut) {
        double x0, y0; // vector origin
        double dx, dy, // vector components
                l; // vector length
        int cIndex;
        //String[] colors = new String[256];
        String[] rgb = new String[3];

        // draw vectors
        for (int i = 0; i < d.length; ++i) {
            ip.setColor(new Color(255, 0, 0));
            x0 = d[i][0];
            y0 = d[i][1];
            dx = d[i][2] * vs;
            dy = d[i][3] * vs;

            l = magArr[i % dim[0]][(int) Math.floor(i / dim[0])];
            //l = Math.sqrt((dx * dx + dy * dy));
//            if (autoScale) {
//                cIndex = (int) (l * 255 / max);
//            } else {
//                cIndex = (int) (l * 255 / colorMax);
//            }
            cIndex = (int) (l * 255 / colorMax);
            if (cIndex < 0) {
                cIndex = 0;
            }
            if (cIndex > 255) {
                cIndex = 255;
            }
            rgb = lut[cIndex].split(" ");
            ip.setColor(new Color(Integer.parseInt(rgb[0]), Integer.parseInt(rgb[1]), Integer.parseInt(rgb[2])));
            Arrow a = new ij.gui.Arrow(x0, y0, x0 + dx, y0 + dy);
            if (l * vs / 3.5 > 4) {
                a.setHeadSize(l * vs / 4);
            } else {
                a.setHeadSize(4);
            }

            if(dx!=0 || dy!=0) ip.draw(a);
        }
    }

    private static String[] getLUT(String lut) throws Exception {

        String[] lines;
        String[] RGB = new String[256];
        String path = IJ.getDirectory("luts") + lut + ".lut";
        int j;

        try {
            String file = IJ.openAsString(path);
        } catch (Exception e) {
            IJ.error("error opening lut file");
        }
        String file = IJ.openAsString(path);
        lines = file.split("\n");
        if (lines.length < 255 || lines.length > 257) {
            throw new Exception("unsupported LUT format");
        }
        String firstStr = lines[0].substring(0, 1);
        try {
            int f = Integer.parseInt(firstStr);
            j = 0;
        } catch (Exception e) {
            j = 1;
        }
        for (int i = j; i < 256; i++) {
            RGB[i] = lines[i];
        }
        return RGB;

    }

    private static double[][] calcMag(double[][][] vectorArr) {

        double[][] Arr = new double[vectorArr.length][vectorArr[0].length];

        for (int j = 0; j < vectorArr.length; j++) {
            for (int i = 0; i < vectorArr[0].length; i++) {

                Arr[j][i] = Math.sqrt(vectorArr[j][i][2] * vectorArr[j][i][2] + vectorArr[j][i][3] * vectorArr[j][i][3]);

            }
        }

        return Arr;
    }

    protected static double[][][] convert2DPivTo3D(double[][] _piv, int nx, int ny) {

        //check the _piv array is fixing x first or y first
        boolean fx;
        if (_piv[0][0] == _piv[1][0]) // x coordinates is fixed first
        {
            fx = true;
        } else {
            fx = false;
        }

        double[][][] newPIV = new double[nx][ny][_piv[0].length];
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                for (int k = 0; k < _piv[0].length; k++) {
                    if (fx) {
                        newPIV[i][j][k] = _piv[ny * i + j][k];
                    } else {
                        newPIV[i][j][k] = _piv[nx * j + i][k];
                    }
                }
            }
        }

        return newPIV;

    }

    protected static double findMax2DArray(double[][] M) {

        double m = M[0][0];

        for (int i = 0; i < M.length; i++) {
            for (int j = 0; j < M[0].length; j++) {
                if (M[i][j] > m) {
                    m = M[i][j];
                }
            }
        }
        return m;
    }

    static protected float[][] doubleArrayToFloat(double[][] dArr) {

        float[][] fArr = new float[dArr.length][dArr[0].length];

        for (int i = 0; i < fArr.length; i++) {
            for (int j = 0; j < fArr[0].length; j++) {
                fArr[i][j] = (float) dArr[i][j];
            }
        }

        return fArr;
    }

    private static void splitLUTtoRGB(String[] _lut){

        String[] rgb = new String[3];
        lutR = new byte[_lut.length];
        lutG = new byte[_lut.length];
        lutB = new byte[_lut.length];

        for(int i=0;i<_lut.length;i++){
            rgb =  _lut[i].split(" ");
            lutR[i] = (byte)Integer.parseInt(rgb[0]);
            lutG[i] = (byte)Integer.parseInt(rgb[1]);
            lutB[i] = (byte)Integer.parseInt(rgb[2]);

        }

        return;
    }

   protected static void loadLut(String lut){

        if (!lut.equals("S_Pet")) {
            try {
                colors = getLUT(lut);
            } catch (Exception e) {
                IJ.error("error reading LUT file");
            }
        } else {
            colors = spet;
        }
        splitLUTtoRGB(colors);
        cm = new IndexColorModel(8, 256, lutR, lutG, lutB);
    }


    protected void calculateMag(double[][] _d){

        double[][][] matrix;
        //convert the 2D ([#data points][x,y,dx,dy]) PIV array to 3D ([#data points in width][#data points in height][x,y,dx,dy])
        matrix = convert2DPivTo3D(_d, dimensions[0], dimensions[1]);
        magArr = calcMag(matrix);
    }


    protected static double[][] load2DArrayFromFile(String path) throws Exception{

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
            //IJ.error(e.getMessage());
            throw new Exception("Unsupported file format.");
        }

        matrix = new double[row][5];
        Iterator<String[]> iter = aL.iterator();
        int counter = 0;

        /* go over all elements of aL*/
        while (iter.hasNext()) {
            cell = (String[]) iter.next();     // cell correspond to one element(string[4]) of aL, that is: one entry of data containing x,y,ux,uy, L
            //check if we have the fifth column which correspond to the vector magnitude in our data.
            if (cell.length>4){
                for (int i = 0; i < 5; i++) {
                    //matrix[counter][i] = Double.parseDouble(cell[i]);
                    matrix[counter][i] = nf.parse(cell[i]).doubleValue();
                }
                counter++;
            }else if (cell.length==4){
                for (int i = 0; i < 4; i++) {
                    //matrix[counter][i] = Double.parseDouble(cell[i]);
                    matrix[counter][i] = (Double)nf.parse(cell[i]).doubleValue();
                }
                matrix[counter][4] = Math.sqrt(matrix[counter][2]*matrix[counter][2]+matrix[counter][3]*matrix[counter][3]);
                counter++;
            }else{
                throw new Exception("The file must have at least first 4 column: x,y,dx,dy separated by space or tab");
            }
        }
        return matrix;
    }

   /**
     *extract one element from a 2D array [nData points][nElements]
     *into another 2D data array with this specified element indexed by
     * [nData points in x][nData points in y]
     */
   protected static double[][] get2DElement(double[][] _piv, int nx, int ny, int nEle) {

        //check the _piv array is fixing x first or y first
        boolean fx;
        if (_piv[0][0] == _piv[1][0]) // x coordinates is fixed first
        {
            fx = true;
        } else {
            fx = false;
        }

        double[][] new2D = new double[nx][ny];
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                if (fx) {
                    new2D[i][j] = _piv[ny * i + j][nEle];
                } else {
                    new2D[i][j] = _piv[nx * j + i][nEle];
                }
            }
        }

        return new2D;

    }

}
