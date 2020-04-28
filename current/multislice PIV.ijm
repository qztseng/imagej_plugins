id0 = getImageID();
slices = nSlices;

Dialog.create("Multi-slice PIV");
Dialog.addNumber("Interrogation window 1", 128);
Dialog.addNumber("search window 1", 256);
Dialog.addNumber("vector spacing 1", 64);
Dialog.addNumber("Interrogation window 2", 64);
Dialog.addNumber("search window 2", 128);
Dialog.addNumber("vector spacing 2", 32);
Dialog.addNumber("Interrogation window 3", 48);
Dialog.addNumber("search window 3", 96);
Dialog.addNumber("vector spacing 3", 24);
Dialog.addNumber("correlation threshold", 0.8);
Dialog.addCheckbox("Always match with 1st slice?", false);
Dialog.show;

piv1 = Dialog.getNumber();
sw1 = Dialog.getNumber();
vs1 = Dialog.getNumber();
piv2 = Dialog.getNumber();
sw2 = Dialog.getNumber();
vs2 = Dialog.getNumber();
piv3 = Dialog.getNumber();
sw3 = Dialog.getNumber();
vs3 = Dialog.getNumber();
corr = Dialog.getNumber();
first = Dialog.getCheckbox();

path = getDirectory("select a folder to save PIV results");

setBatchMode(true);
for(s=1;s<slices;s++){
	selectImage(id0);
	pad = floor(log(slices)/log(10))+1;
	spad = IJ.pad(s,pad);
	if(first){
		setSlice(1);
		run("Duplicate...", "title=FIRST range=1-1");
		selectImage(id0);
		setSlice(s+1);
		run("Duplicate...", "title="+(s+1)+" range="+(s+1)+"");
		run("Concatenate...", "  image1=FIRST image2="+(s+1)+" image3=[-- None --]");
		rename("seq_"+spad);
	}else{
		run("Duplicate...", "title=[seq_"+spad+"] duplicate range="+s+"-"+(s+1)+"");
	}
	run("iterative PIV(Advanced)...", "  piv1="+piv1+" sw1="+sw1+" vs1="+vs1+" piv2="+piv2+" sw2="+sw2+" vs2="+vs2+" piv3="+piv3+" sw3="+sw3+" vs3="+vs3+" correlation="+corr+" batch path=["+path+"]");
	nI = nImages;
	//print(nI);
	for(n=1;n<nI;n++){
		selectImage(1);
		//print(n);
		//print(getTitle());
		if(getImageID!=id0) close();
	}
}
		
