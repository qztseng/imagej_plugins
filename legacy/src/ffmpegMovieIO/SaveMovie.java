package ffmpegMovieIO;

import ij.*;
import ij.process.*;
import ij.plugin.filter.*;
import ij.plugin.*;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.gui.GenericDialog;
import ij.io.SaveDialog;
import java.awt.image.BufferedImage;
import static org.bytedeco.javacpp.opencv_core.*;
import static org.bytedeco.javacpp.avcodec.*;
import static org.bytedeco.javacpp.avutil.*;
import org.bytedeco.javacv.FFmpegFrameRecorder;
import org.bytedeco.javacv.FrameRecorder;

/**This is a video output plugin for ImageJ. The video encoding function is using the opencv, and ffmpeg libraries through javacv and javacpp. 
 * Currently only MPEG4 (MPEG-4 part 2),FLV (Flash Video (FLV) / Sorenson Spark / Sorenson H.263), WMV(Windows Media Video 7), 
 * and H.264(libx264 H.264 / AVC / MPEG-4 AVC / MPEG-4 part 10) codecs as well as uncompressed rawvideo are implemented.
 *   
 * First version 2012/6/12
 * Second version 2012/6/17  add H.264 support
 * Third version 2012/7/17 add option to reduce memory usage
 * v4	2015/02/06 update to the latest javacv 0.10
 * @author Qingzong TSENG, qztseng "at" gmail com. 
 */
public class SaveMovie implements PlugInFilter {
  
    ImagePlus imp1;
    static int container = -1;
    String[] containers = {".avi", ".mov",".mp4"};
    static int codec = -1;
    String[] codecs = {"MPEG4", "RAW", "FLV", "WMV7", "H.264"};
    static int quality = -1;
    String[] qualities = {"custom", "low", "normal", "high", "excellent"};
    int[] codecsID = {AV_CODEC_ID_MPEG4, AV_CODEC_ID_RAWVIDEO, AV_CODEC_ID_FLV1, AV_CODEC_ID_WMV1, AV_CODEC_ID_H264};
    String dir, file, title;
    double fps = 24;
    int br = 400;
    int width, height;
    String ct;
    boolean speed = false;
    boolean zSerie = false, tSerie = false;
    int[] dim;
    
    @Override
    public int setup(String arg, ImagePlus imp) {
        
        if(arg.equals("about")){showAbout();return DONE;}
        if(IJ.getVersion().compareTo("1.46d")<0){
            IJ.log("ImageJ 1.46d or higher version required");
            return DONE;
        }
        
        try {
            imp1=imp;
            title = imp.getTitle();
            dim = imp.getDimensions();
            if(dim[3]!=1){
                if(dim[4]!=1){
                    IJ.showMessage("multi-dimensional stack not supported\n (Reduce dimension to either single Z-slice or single time-point)");
                    return DONE;
                }else zSerie = true;
            }else if(dim[3]==1){
                if(dim[4]!=1){
                    tSerie = true;
                }else{
//                    IJ.showMessage("multi-slice or multi-frame stack required");
//                    return DONE;
                }
            }    

            //if(imp.getType()!= ImagePlus.COLOR_RGB) toRGB = true;
            //if(imp.isComposite() || imp.getOverlay()!=null) toRGB = true;
        } catch (NullPointerException ex) {
            IJ.showMessage("Image stack required");
            return DONE;
        }
        
        return DOES_ALL|STACK_REQUIRED ;
    }
    
    public void run (String arg){
        
        width = imp1.getWidth();
        height = imp1.getHeight();
        ImageProcessor  ipOri;
        IplImage image;
        int frames;
        BufferedImage bi;
        
        if(arg.equals("about")){showAbout();return;}
        if(!showDialog()) return;
        
        if(zSerie) frames = imp1.getNSlices();
        else if (tSerie) frames = imp1.getNFrames();
        else frames = imp1.getNChannels();
        
        //configure the FFmpeg recorder
        FFmpegFrameRecorder recorder = new FFmpegFrameRecorder(dir+file, width, height);
        recorder.setVideoCodec(codecsID[codec]);
        recorder.setFormat(ct);
        recorder.setPixelFormat(AV_PIX_FMT_YUV420P);
        recorder.setFrameRate(fps);
        recorder.setVideoBitrate(br);
        
        
        if(speed){
            try {
                recorder.start();
                
                for (int i = 1; i <= frames; i++) {
                    IJ.showStatus("saving movie...");
                    IJ.showProgress(i+1, frames);
                    if(zSerie)          {imp1.setZ(i);}
                    else if (tSerie)    {imp1.setT(i);}
                    else                {imp1.setC(i);}

                    bi = (imp1.flatten()).getBufferedImage();
                    image = IplImage.createFrom(bi);
                    recorder.record(image, AV_PIX_FMT_BGRA);
                    image.release();
                    bi.flush();
                }
                recorder.stop();
            } catch (FrameRecorder.Exception e) {
                IJ.log("error:"+e);
            }
            
            
        }else {
            ImageStack stk = getFlattenRGB(imp1);
            try {
                recorder.start();
                for (int i = 1; i <= stk.getSize(); i++) {
                    IJ.showStatus("saving movie...");
                    IJ.showProgress(i+1, frames);
                    ipOri = stk.getProcessor(i);
                    bi = ipOri.getBufferedImage();
                    image = IplImage.createFrom(bi);
                    recorder.record(image, AV_PIX_FMT_BGRA);
                    image.release();
                    bi.flush();
                }
                recorder.stop();
            } catch (FrameRecorder.Exception e) {
                IJ.log("error:"+e);
            }
            
         
        }
    }
    
    
    private ImageStack getFlattenRGB(ImagePlus imp){
        
        ImagePlus imp2;
        Overlay OL = imp.getOverlay();
        ImageStack stk;
        int ss = imp.getNSlices();
        int ff = imp.getNFrames();
        
        
        if(imp.isComposite()){
            imp2 = imp.createHyperStack("", 1, ss, ff, 24);
            new RGBStackConverter().convertHyperstack(imp,imp2);
            if(OL != null){
                imp2.setOverlay(OL);
                IJ.run(imp2, "Flatten", "stack");
            }
            //imp2.show();
            stk = imp2.getStack();
            imp2.close();
        }else if (imp.getType()== ImagePlus.COLOR_RGB && OL == null){
            stk = imp.getStack();
        }else{
            imp2 = imp.duplicate();
            new StackConverter(imp2).convertToRGB();
            if(OL != null){
                imp2.setOverlay(OL);
                IJ.run(imp2, "Flatten", "stack");
            }
            stk = imp2.getStack();
            imp2.close();
        }
        
        
//        if (dim[3] != 1) {
//            IJ.run(imp2, "RGB Color", "slics");IJ.log("slice");
//            imp2.show();
//        } else if (dim[4] != 1) {
//            IJ.run(imp2, "RGB Color", "frames");IJ.log("frames");
//        } else {
//            IJ.run(imp2, "RGB Color", "");IJ.log("??");
//        }
//        if (imp2.getOverlay() != null) {
//            IJ.run(imp2, "Flatten", "stack");IJ.log("flatten");
//        }
        //ImageStack stk = imp2.getStack();
        //imp2.close();
       
//        }else if (imp.getOverlay()!=null){
//            IJ.log("Flatten2");
//            if(imp.getType()!= ImagePlus.COLOR_RGB){            //overlay flatten only works on RGB stack
//                new StackConverter(imp).convertToRGB();
//            }
//            Overlay overlay = imp.getOverlay();
//            for (int i=1; i<=stk.getSize(); i++) {
//                ImageProcessor OLip = stk.getProcessor(i);
//                Roi[] rois = overlay.toArray();
//                for (int j=0; j<rois.length; j++) {
//                    Roi roi = rois[j];
//                    int position = roi.getPosition();
//                    if (position==0 || position==i)
//                        OLip.drawRoi(roi);
//                }
//            }   
//           
//            imp.setOverlay(null);
//        
//            imp2=imp.flatten();
//            IJ.log("nSlices:"+imp2.getNSlices());
//            IJ.log("nFrames:"+imp2.getNFrames());
//            if(imp2.getType()== ImagePlus.COLOR_RGB) IJ.log("COLOR_RGB now");
//            stk = imp2.getStack();
//            imp2.close();
//        }else if (imp.getType()!= ImagePlus.COLOR_RGB) {
//            IJ.log("Flatten3");
//            new StackConverter(imp).convertToRGB();
//            stk = imp.getStack();
//        }else{
//            IJ.log("Flatten4");
//        }    
        return stk;
    }
    
    public boolean showDialog() {
    	
        String defaultItem;
        
        GenericDialog gd = new GenericDialog("saveMovie", IJ.getInstance());
        
        
        fps = Animator.getFrameRate();
        gd.addNumericField("Frame rate", fps, 1);
        if (container == -1) {
            defaultItem = containers[0];
        } else {
            defaultItem = containers[container];
        }
        gd.addChoice("container format", containers, defaultItem);
       // gd.addStringField("container format", "avi");
        
        if (codec == -1) {
            defaultItem = codecs[0];
        } else {
            defaultItem = codecs[codec];
        }
        gd.addChoice("using codec", codecs, defaultItem);
        
        if (quality == -1) {
            defaultItem = qualities[2];
        } else {
            defaultItem = qualities[quality];
        }
        gd.addChoice("video quality", qualities, defaultItem);
        gd.addNumericField("custom bitrate(kb/s)", br, 0);
        gd.addCheckbox("Use less memory (slower)", false);
        gd.showDialog();
        
        if (gd.wasCanceled()) {
            
            return false;
        }
        
        fps = gd.getNextNumber();
        container = gd.getNextChoiceIndex();
        //ct = gd.getNextString();
        codec = gd.getNextChoiceIndex();
        
        ct = containers[container].substring(1);
        if (codecs[codec].equals("FLV")) {
            ct = "flv";
        }
        if (codecs[codec].equals("RAW")) {
            ct = "avi";
        }
        if (codecs[codec].equals("WMV7") && (ct.equals("mov") || ct.equals("mp4"))) {
            ct = "wmv";
        }
        if (codecs[codec].equals("H.264") ) {
            //IJ.showMessage("warning", "H.264 codec requires installing FFmpeg with libx264 support");
            ct = "mov";
        }
        
        quality = gd.getNextChoiceIndex();
        
        if(quality==0)br = (int)gd.getNextNumber() * 1000;
        else br = estimateBitrate(width, height, fps, quality);
        
        speed = gd.getNextBoolean();
        
        //SaveDialog sd = new SaveDialog("Save Video", IJ.getDirectory("home"), imp.getTitle(), containers[container]);
        SaveDialog sd = new SaveDialog("Save Video", IJ.getDirectory("home"), title, "."+ct);
        if (sd.getDirectory() == null || sd.getFileName() == null) {
            return false;
        }
        dir = sd.getDirectory();
        file = sd.getFileName();
        
        return true;
    
    }
    
    protected int estimateBitrate (int width, int height, double FPS, int mFactor){
        
        int bitrate = (int)Math.round(width*height*FPS*mFactor*0.1);   
        
        return bitrate;
    }
    
     private void showAbout() {
		IJ.showMessage("save as Movie", "This is a video output plugin for ImageJ. The video encoding function is using \n"
                        + "the opencv and ffmpeg libraries through javacv and javacpp. Currently only MPEG4 \n"
                        + "(MPEG-4 part 2), FLV (Flash Video (FLV) / Sorenson Spark / Sorenson H.263).\n"
                        + "WMV (Windows Media Video 7), and H.264 (libx264 H.264 / AVC / MPEG-4 AVC / MPEG-4 part 10)\n"
                        + "codecs as well as uncompressed rawvideo are implemented.\n"
                        + "Please send comments to qztseng at gmail dot com"
                        
			);
	}

    @Override
    public void run(ImageProcessor ip) {
        
        run("");
    }
    
  
    
}
