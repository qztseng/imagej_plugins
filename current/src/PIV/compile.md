  javac -source 1.8 -target 1.8 -d . -cp /Applications/ImageJ.app/Contents/Java/ij.jar:/Applications/ImageJ.app/plugins/jars/*:. *.java 

  jar -cvf ../../PIV_.jar ./plugins.config ./PIV/*.class ./org

  jar -uvf ../../PIV_.jar PIV/*.class
