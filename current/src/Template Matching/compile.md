javac -source 1.8 -target 1.8 -d . -cp $PATH_OF_IJ.JAR/ij.jar:$PATH_OF_JAVACV_JARS/*:. *.java
jar -cvf Template_Matching.jar plugins.config ./TemplateMatching/

$PATH_OF_IJ.JAR and $PATH_OF_JAVACV_JARS should be the folder where you have the ij.jar and the opencv/javacv library jars. 
For example: /Applications/ImageJ.app/Contents/Java/ij.jar and /Applications/ImageJ.app/plugins/jars
The plugins.config is a special configuration file for ImageJ plugin, you can get it by simply extract the Template_Matching.jar.

If you are running it under windows system, you have to replace the colon (:) with semi-colon(;)
