package LPMatching;

import ij.*;
import ij.process.*;
import ij.process.FHT.*;
import ij.gui.*;
import java.awt.Color;
import ij.plugin.*;
import ij.plugin.filter.MaximumFinder;
import ij.measure.ResultsTable;
import java.awt.image.BufferedImage;
import java.nio.FloatBuffer;
import java.util.List;
import java.util.*;
//import java.lang.RuntimeException;
import org.bytedeco.javacv.*;
import org.bytedeco.javacpp.*;
import org.bytedeco.javacv.Java2DFrameUtils;
import static org.bytedeco.javacpp.opencv_core.*;
import static org.bytedeco.javacpp.opencv_imgproc.*;


/*
 * This plugin implement the template match method from OpenCV library. The interface between OpenCV and Java is 
 * using Samuel Audet's JavaCV code from: http://code.google.com/p/javacv/
 * The detailed algorithm of each matching method can be found at opencv's documentation page:
 * http://opencv.willowgarage.com/documentation/object_detection.html?highlight=matchtemplate#cvMatchTemplate
 * It supports 8 bit and 16 bit grayscale only. 
 * 
 * 2013/05/14 Add support for 32bit grayscale
 * 2015/02/06 Update to latest javacv v0.10
 * 2019/03/17 Update to work with javacv 1.4.4
 * 
 * By TSENG Qingzong (qztseng /at/ gmail dot com) 
 */
public class logpolar_match implements PlugIn {

    private static String title1 = "";
    private static String title2 = "";
    private static String[] titles;
    private static boolean showRT = false;
    private static boolean log = true;
    private double magnitude = 56;
    int[] dxdy = new int[2];
    private static ImagePlus img1, img2;
    private BufferedImage bi, bi2;
    private int[] wList;
	private ByteProcessor p1,p2;
	private FloatProcessor m1,m2;

	
    @Override
    public void run(String arg) {

        if (arg.equals("about")) {
            showAbout();
            return;
        }

        /*
         ** check how many existing images
         */

        wList = WindowManager.getIDList();
        if (wList == null) {
            IJ.noImage();
            return;
        }
        int wCount = WindowManager.getImageCount();
        if (wCount < 2) {
            IJ.error("we need two images to do template matching");
            return;
        }

        /*
         ** show Dialog and get user input parameters
         */


        if (!getParams()) {
            return;
        }

        /*
         ** check parameter compatibility
         */

        if (img1.getBitDepth() != img2.getBitDepth()) {
            IJ.error("Images need to have the same type (bit depth)");
            return;
        }

        /*
         ** start matching
         */
		Loader.load(opencv_core.class);

        long startTime = System.currentTimeMillis();

        Point2d pt = doMatchLP(img1, img2, magnitude);

        long elapsedTime = System.currentTimeMillis() - startTime;

        double scale = Math.exp( pt.x() / magnitude);
		double rotate = pt.y()*360/img1.getWidth();

		IJ.log("pt.x:"+pt.x()+", pt.y:"+pt.y()+", scale:"+scale+", rotate:"+rotate);
		
		new ImagePlus("p1", p1).show();
		new ImagePlus("p2", p2).show();
		new ImagePlus("m1", m1).show();
		new ImagePlus("m2", m2).show();
		//new ImagePlus("FC", phaseC).show();
		
        IJ.showStatus("Matching done in " + elapsedTime + "msec");
        //IJ.log("Matching done in "+elapsedTime+"msec");

        return;
    }

    private boolean getParams() {

        titles = new String[wList.length];
        for (int i = 0; i < wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
            if (imp != null) {
                titles[i] = imp.getTitle();
            } else {
                titles[i] = "";
            }
        }

        GenericDialog gd = new GenericDialog("MatchTemplate", IJ.getInstance());
        String defaultItem;

        if (title1.equals("")) {
            defaultItem = titles[0];
        } else {
            defaultItem = title1;
        }
        gd.addChoice("Image:", titles, defaultItem);
        if (title2.equals("")) {
            defaultItem = titles[1];
        } else {
            defaultItem = title2;
        }
        gd.addChoice("Template:", titles, defaultItem);
        gd.addNumericField("Magnitude for LogPolar transform", magnitude, 0);
		gd.addCheckbox("Output in resultTable ?", showRT);
        gd.addCheckbox("log result?", log);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return false;
        }
        int index1 = gd.getNextChoiceIndex();
        title1 = titles[index1];
        int index2 = gd.getNextChoiceIndex();
        title2 = titles[index2];
        img1 = WindowManager.getImage(wList[index1]);
        img2 = WindowManager.getImage(wList[index2]);
        magnitude = gd.getNextNumber();
		showRT = gd.getNextBoolean();
        log = gd.getNextBoolean();

        return true;
    }

    public Point2d doMatchLP(ImagePlus src, ImagePlus tpl, double magnitude) {
        
        return doMatchLP(src.getProcessor(), tpl.getProcessor(), magnitude);

    }

    public Point2d doMatchLP(ImageProcessor src, ImageProcessor tpl, double magnitude) {

        BufferedImage bi = null, bi2 = null;
        FloatProcessor resultFp = null;
        int srcW = src.getWidth();
        int srcH = src.getHeight();
        int tplW = tpl.getWidth();
        int tplH = tpl.getHeight();
        IplImage ipl_Src, ipl_Tpl, ipl_SrcP, ipl_TplP;
		//Mat SrcP, TplP, SrcP64=null, TplP64=null;
        Mat Src = new Mat(), Tpl = new Mat(), SrcP64 = new Mat(), TplP64 = new Mat(), SrcM = new Mat(), TplM = new Mat();
		 
        switch (src.getBitDepth()) {

/*             case 32:
                //convert source imageProcessor into iplImage
                CvMat srcMat = CvMat.create(srcH, srcW, CV_32FC1);
                double[] dArr1 = float2DtoDouble1DArray(src.getFloatArray(), srcW, srcH);
                srcMat.put(0, dArr1, 0, dArr1.length);
                iplSrc = srcMat.asIplImage();
                //iplSrc = temp.clone();
                
                //convert template imageProcessor into iplImage
                CvMat tplMat = CvMat.create(tplH, tplW, CV_32FC1);
                double[] dArr2 = float2DtoDouble1DArray(tpl.getFloatArray(), tplW, tplH);
                tplMat.put(0, dArr2, 0, dArr2.length);
                iplTpl = tplMat.asIplImage();
                //iplTpl = temp2.clone();
                
                break;
            case 16:
                //since cvMatchTemplate don't accept 16bit image, we have to convert it to 32bit
                iplSrc = cvCreateImage(cvSize(srcW, srcH), IPL_DEPTH_32F, 1);
                bi = ((ShortProcessor) src).get16BitBufferedImage();
                //temp = IplImage.createFrom(bi);
                temp = Java2DFrameUtils.toIplImage(bi);
                cvConvertScale(temp, iplSrc, 1 / 65535.0, 0);

                iplTpl = cvCreateImage(cvSize(tplW, tplH), IPL_DEPTH_32F, 1);
                bi2 = ((ShortProcessor) tpl).get16BitBufferedImage();
                //temp2 = IplImage.createFrom(bi2);
                temp2 = Java2DFrameUtils.toIplImage(bi2);
                cvConvertScale(temp2, iplTpl, 1 / 65535.0, 0);

                temp.release();
                temp2.release();

                break; */
            case 8:

                bi = src.getBufferedImage();
				ipl_Src = Java2DFrameUtils.toIplImage(bi);
				Src = Java2DFrameUtils.toMat(src.getBufferedImage());
                bi2 = tpl.getBufferedImage();
                ipl_Tpl = Java2DFrameUtils.toIplImage(bi2);
				Tpl = Java2DFrameUtils.toMat(tpl.getBufferedImage());

                break;
            default:
                //ipl_Src = null;
				//ipl_Tpl = null;
				IJ.error("Unsupported image type");
                throw new RuntimeException("Unsupported Image Bit Depth:" + src.getBitDepth());
				//break;
        }

		ipl_SrcP = Java2DFrameUtils.toIplImage(Mat.zeros(srcH, srcW, CV_8UC1).asMat());
		ipl_TplP = Java2DFrameUtils.toIplImage(Mat.zeros(tplH, tplW, CV_8UC1).asMat());
		
		FFTMagnitude(Src, SrcM);
		FFTMagnitude(Tpl, TplM);
		
		
		//FloatBuffer fb = SrcM.<FloatBuffer>createBuffer();
		FloatBuffer fb = SrcM.createBuffer();
/*         float[] f = new float[SrcM.cols() * SrcM.rows()];
        fb.get(f, 0, f.length);
        m1 = new FloatProcessor(SrcM.cols(), TplM.rows(), f, null); */
		m1 = new ByteProcessor(fb);
		
		fb = TplM.<FloatBuffer>createBuffer();
		f = new float[TplM.cols() * TplM.rows()];
		fb.get(f, 0, f.length);
        m2 = new FloatProcessor(TplM.cols(), TplM.rows(), f, null);

		
		//ipl_SrcP = SrcP;
		//ipl_TplP = TplP;
		
		cvLogPolar(ipl_Src, ipl_SrcP, cvPoint2D32f(srcW/2, srcH/2), magnitude);
		p1 = new ByteProcessor(Java2DFrameUtils.toBufferedImage(ipl_SrcP));
			
		cvLogPolar(ipl_Tpl, ipl_TplP, cvPoint2D32f(tplW/2, srcH/2), magnitude);
		p2 = new ByteProcessor(Java2DFrameUtils.toBufferedImage(ipl_TplP));
		
		Java2DFrameUtils.toMat(ipl_SrcP).convertTo(SrcP64, CV_64F);
		Java2DFrameUtils.toMat(ipl_TplP).convertTo(TplP64, CV_64F);
		
        Point2d pt = phaseCorrelate(SrcP64, TplP64);

/*         cvMatchTemplate(iplSrc, iplTpl, res, method);
        FloatBuffer fb = res.getFloatBuffer();
        float[] f = new float[res.width() * res.height()];
        fb.get(f, 0, f.length);
        resultFp = new FloatProcessor(res.width(), res.height(), f, null);
        cvReleaseImage(res);
        
        switch (src.getBitDepth()) {
            case 32:
                iplSrc.release();
                iplTpl.release();
                break;
            case 16:
                cvReleaseImage(iplSrc);
                cvReleaseImage(iplTpl);
                break;
            case 8:
                iplSrc.release();
                iplTpl.release();
                break;
            default:
                break;    
        } */

        return pt;
    }
	
	public void FFTMagnitude(Mat Src, Mat magI){
		
		/*
		Code copied from https://docs.opencv.org/3.4/d8/d01/tutorial_discrete_fourier_transform.html
		*/
		int M = getOptimalDFTSize(Src.rows());
		int N = getOptimalDFTSize(Src.cols());
		
		Mat padded = new Mat();
		copyMakeBorder(Src, padded, 0, M - Src.rows(), 0, N - Src.cols(), BORDER_CONSTANT, Scalar.all(0));
		
		magI = Src;
		
/* 		MatVector planes = new MatVector(2);
        padded.convertTo(padded, CV_32F);
        planes.put(0, padded);
        planes.put(1, Mat.zeros(padded.size(), CV_32F).asMat());
        
		Mat complexI = new Mat();
        merge(planes, complexI);         // Add to the expanded another plane with zeros
        dft(complexI, complexI);         // this way the result may fit in the source matrix
        // compute the magnitude and switch to logarithmic scale
        // => log(1 + sqrt(Re(DFT(I))^2 + Im(DFT(I))^2))
        split(complexI, planes);                               // planes.get(0) = Re(DFT(I)
                                                                    // planes.get(1) = Im(DFT(I))
        magnitude(planes.get(0), planes.get(1), planes.get(0));// planes.get(0) = magnitude
        magI = planes.get(0); */
		
		return;
	}

    public void showAbout() {
        IJ.showMessage("cvMatch Template", "This plugin implements the tempalte matching function from\n"
                + "the OpenCV library. It will try to find an object (template)\n"
                + "within a given image (image).Six different matching algorithms\n"
                + "(methods)could be used. The matching result could be printed\n"
                + "in the log window or the result table. \n"
                + "You can also decide to display the correlation map. The coordinates\n"
                + "of the maximum (or minimum for the square difference methods)\n"
                + "correspond to the best match.\n"
                + "By checking the multimatch option, not only the best match will\n"
                + "be shown, but also all the similar pattern above the defined\n"
                + "threshold will be shown (Find maximum function on the correlation map)\n"
                + "More details on \nhttps://sites.google.com/site/qingzongtseng/template-matching-ij-plugin");
    }
}
