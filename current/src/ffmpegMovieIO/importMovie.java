package ffmpegMovieIO;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.io.OpenDialog;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import java.awt.image.BufferedImage;
import static org.bytedeco.javacpp.opencv_core.*;
import static org.bytedeco.javacpp.opencv_imgproc.*;
import static org.bytedeco.javacpp.avcodec.*;
import static org.bytedeco.javacpp.avutil.*;
import org.bytedeco.javacv.*;
import org.bytedeco.javacpp.*;
import org.bytedeco.javacv.FFmpegFrameGrabber;
import org.bytedeco.javacv.FrameGrabber;

/**This is a video output plugin for ImageJ. The video encoding function is using the opencv, and ffmpeg libraries through javacv and javacpp. 
 * Currently only MPEG4 (MPEG-4 part 2),FLV (Flash Video (FLV) / Sorenson Spark / Sorenson H.263), WMV(Windows Media Video 7), 
 * and H.264(libx264 H.264 / AVC / MPEG-4 AVC / MPEG-4 part 10) codecs as well as uncompressed rawvideo are implemented.
 *   
 * First version 2012/6/12
 * Second version 2012/6/17  add H.264 support
 * Third version 2012/7/17 add option to reduce memory usage
 * v4	2015/02/06 update to the latest javacv 0.10
 * v5   2019/03/18 update to work with javacv 1.4.4
 * @author Qingzong TSENG, qztseng "at" gmail com. 
 */

public class importMovie
  implements PlugIn
{
  String dir;
  String file;
  String title;
  int width;
  int height;
  int pix_fmt;
  int nFrame;
  int startFrame;
  int endFrame;
  int counter;
  double frameRate;
  boolean virtual;
  boolean convertGray;

  public void run(String arg)
  {
    if (arg.equals("about")) { showAbout(); return; }
    if (!showDialog()) return;

    FrameGrabber fG = new FFmpegFrameGrabber(this.dir + this.file);
	OpenCVFrameConverter.ToIplImage converter = new OpenCVFrameConverter.ToIplImage();
    IplImage cFrame2;
	BufferedImage bi, bi2;

    this.counter = 1;
    try
    {
      fG.start();
	  Frame cFrame = null;
	  
      this.width = fG.getImageWidth();
      this.height = fG.getImageHeight();
      this.pix_fmt = fG.getPixelFormat();
      this.nFrame = fG.getLengthInFrames();
      this.frameRate = fG.getFrameRate();

      if (!sDialog2()) return;

      ImageStack is = new ImageStack(this.width, this.height);

      while (this.counter <= this.endFrame) {
        this.counter += 1;
        try {
          cFrame = fG.grab();
          if (cFrame == null) {
            IJ.log("counter:" + this.counter);
            break;
          }
          if (this.counter > this.startFrame)
          {
            bi = new Java2DFrameConverter().convert(cFrame);

            if (bi.getType() == 10) {
              ByteProcessor bp = new ByteProcessor(bi);

              is.addSlice(bp);
            }
            else if (!this.convertGray) {
              ColorProcessor cp = new ColorProcessor(bi);

              is.addSlice(cp);
            } else {
              cFrame2 = IplImage.create(cFrame.imageWidth, cFrame.imageHeight, cFrame.imageDepth, 1);
              cvCvtColor(converter.convert(cFrame), cFrame2, 6);
              bi2 = Java2DFrameUtils.toBufferedImage(cFrame2);
              ByteProcessor bp = new ByteProcessor(bi2);
              is.addSlice(bp);
            }
          }
        }
        catch (Exception e)
        {
        }
      }
      ImagePlus imp = new ImagePlus(this.title, is);
      imp.show();
      if (this.frameRate != 0.0D) {
        Calibration cal = imp.getCalibration();
        cal.fps = this.frameRate;
        imp.setCalibration(cal);
      }
      fG.release();
    }
    catch (Exception e) {
      IJ.error("Error:" + e);
    }
  }

  private boolean showDialog()
  {
    OpenDialog od = new OpenDialog("Import Movie...", "");
    if (od.getFileName() == null) {
      return false;
    }

    this.dir = od.getDirectory();
    this.file = od.getFileName();
    this.title = this.file.substring(0, this.file.indexOf(".", 0));

    return true;
  }

  private boolean sDialog2()
  {
    GenericDialog gd = new GenericDialog("Import Movie", IJ.getInstance());
    gd.addMessage("Total " + this.nFrame + " frames in the movie.\nSpecify the range of frames you want to import.");
    gd.addNumericField("Starting frame:", 1.0D, 0);
    gd.addNumericField("Endnig frame:", this.nFrame, 0);

    gd.addCheckbox("Convert RGB to 8bit grayscale?", false);
    gd.showDialog();

    this.startFrame = ((int)gd.getNextNumber());
    this.endFrame = ((int)gd.getNextNumber());

    this.convertGray = gd.getNextBoolean();

    return true;
  }

  private void showAbout()
  {
    IJ.showMessage("This is a video import plugin for ImageJ. The video decoding function is using the opencv ,\nand ffmpeg libraries through javacv and javacpp. Theoratically, all the format supported\nby ffmpeg could be imported. Nevertheless, since the setTimeFrame function in ffmpeg is not\nvery stable yet, the support for virtual stack becomes difficult. As as a result, large\nvideos exceeding the memory capacity will have problems. The solution would be import them\npiecewise or use the built-in avi reader by using the avi format. \nPlease send any comment to qztseng /at/ gmail dot com");
  }
}
