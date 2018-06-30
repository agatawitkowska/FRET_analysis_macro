// Macro for quantification of GUV diameter, membrane fluorescence intensity, and FRET
// FRET_analysis_macro.ijm v1.0.0 - based on GUV_membrane_linearization_macro.ijm v1.0.0 doi:10.5281/zenodo.376618
// This macro was developed by Agata Witkowska for the study
// "A convenient protocol for generating giant unilamellar vesicles containing SNARE proteins using electroformation", Witkowska, Jablonski, and Jahn, 2018, Scientific Reports
// https://doi.org/10.1038/s41598-018-27456-4

// This Fiji macro can be used for automatic calculation of FRET donor signal recovery upon acceptor photobleaching in lipid mixing experiments with GUVs (as in Witkowska & Jahn, 2017, doi:10.1016/j.bpj.2017.03.010; and Witkowska, Jablonski, and Jahn, 2018, doi:10.1038/s41598-018-27456-4) in multiple microscopy image stacks (1 GUV on image present). It saves list with membrane fluorescence intensity over time for each image and a Log file with % donor change after photobleaching, file name, and GUV diameter.
// Required plugin: Polar Transformer.

// Generates list of files for analysis
files=newArray(0);
dir=getDirectory("Select the directory"); // here you select a directory where your FRET data files are located
list=getFileList(dir);
for(k=0;k<list.length;k++){
if(endsWith(list[k],".czi")){ // extension of your data files (here Zeiss' ZEN format "*.czi" was used)
        files=Array.concat(files,list[k]);
        }
}
// Open file with Bio-Formats
for(f=0;f<files.length;f++){
       run("Bio-Formats", "open='"+dir+ files[f] + "' color_mode=Default view=Hyperstack stack_order=XYCZT");

// Creates directory for saving files
OrgFile=File.name;
OrgFileName=File.nameWithoutExtension;
File.makeDirectory(File.directory+OrgFileName);
MacroDir=File.directory+OrgFileName+"/";
selectWindow(OrgFile);

// Bleach parameters
PreBleachFrame = getInfo("Experiment|AcquisitionBlock|MultiTrackSetup|TrackSetup|BleachSetup|BleachParameterSet|StartNumber #1"); // retrieve bleach parameters from metadata
PreBleachFrame = parseInt(PreBleachFrame);

// Returns pixel size of original image needed for future diameter calculations
getPixelSize(unit, pixelWidth, pixelHeight);

// Prepare canvas (image) for GUV detection
getDimensions(width, height, channels, slices, frames);
cwidth=width+100;
cheight=height+100;
run("Canvas Size...", "width=cwidth height=cheight position=Center zero"); // width and hight set in this step should exceed analysed image
run("Brightness/Contrast...");
run("Enhance Contrast", "Auto");
selectWindow(OrgFile);
run("Select None");
run("Duplicate...", " ");

// Thresholding and GUV border determination
run("8-bit");
run("Entropy Threshold");
run("Make Binary");
run("Despeckle");
run("Remove Outliers...", "radius=2 threshold=50 which=Dark");
run("Select Bounding Box (guess background color)");
run("Fit Circle");
run("ROI Manager...");
roiManager("Add");
roiManager('select', 0);
roiManager("Rename", "GUV");
close();
selectWindow(OrgFile);
roiManager('select', 0);
waitForUser( "Pause","Check Circle Fit"); // wait for user action until selection is fitted to GUV membrane
roiManager("Update");
selectWindow(OrgFile);

setBatchMode(true);

// Preparation for polar transform
selectWindow(OrgFile);
roiManager('select', 0);
Roi.getBounds(x, y, width, height);
GUV_diameter = width*pixelWidth;
nwidth=width+50;
nheight=height+50;
xC=x+(width/2);
yC=y+(height/2);
run("Specify...", "width=nwidth height=nheight x=xC y=yC oval constrain centered");
roiManager("Add");
roiManager('select', 1);
roiManager("Rename", "GUV_polar");
EditedFileName=OrgFileName+"_ed";
rename(EditedFileName);
run("Duplicate...", "duplicate");
rename(EditedFileName+"_cut");

// Polar transform on a stack
inputId = getImageID(); 
inputTitle = getTitle();
inputName=replace(inputTitle, " ", "_");

Stack.getDimensions(width, height, channels, slices, frames);

for (i=1; i<frames+1; i++){
        Stack.setFrame(i);
        for (j=1; j<channels+1; j++) {
                Stack.setChannel(j);
                run("Polar Transformer", "method=Polar degrees=360 default_center number=360");
                run("Fire");
                selectWindow(inputTitle);
        }
}

run("Images to Stack", "name=PolarStack_"+inputName+" title=[] use");
run("Rotate 90 Degrees Left");
run("Stack to Hyperstack...", "order=xyczt(default) channels="+channels+" slices="+slices+" frames="+frames+" display=Grayscale");
PolarFile =  MacroDir+OrgFileName+"_polar";
PolarFileName = OrgFileName+"_polar";
rename(PolarFileName);

setBatchMode(false);

// Selection of membrane ROI
run("Brightness/Contrast...");
run("Enhance Contrast", "Auto");
selectWindow(PolarFileName);
makeRectangle(0, 100, 360, 20);
roiManager("Add");
roiManager('select', 2);
roiManager("Rename", "membrane");
waitForUser( "Pause","Check Membrane Selection"); // wait for user action until selection is fitted to GUV membrane
roiManager("Update");

// Membrane intensity for channels
selectWindow(PolarFileName);
run("Select None");
run("Duplicate...", "duplicate");
run("Split Channels");
selectWindow("C1-"+PolarFileName+"-1");
roiManager('select', 2);
run("Set Measurements...", "mean standard display redirect=None decimal=3");
roiManager('select', 2);
roiManager("multi-measure measure_all one");
selectWindow("C2-"+PolarFileName+"-1");
roiManager('select', 2);
roiManager("multi-measure measure_all one append");
selectWindow("C1-"+PolarFileName+"-1");
close();
selectWindow("C2-"+PolarFileName+"-1");
close();

// Calculation of the non-bleached channel intensity change upon bleaching
before = 0;
after = 0;
mmin = frames;
mmax = PreBleachFrame+mmin;
for(m=mmin;m<mmax;m++){
before += getResult("Mean(membrane)", m);
}
befmean = before/PreBleachFrame;
for(o=mmax;o<(2*frames);o++){
after += getResult("Mean(membrane)", o);	
}
aftmean = after/(frames-PreBleachFrame);
Ratio = (aftmean/befmean - 1) * 100;

// Log with the channel intensity change (%), file name, and diameter for each analysed file
print(Ratio+"; "+OrgFile+"; "+GUV_diameter);

lastRow = nResults;
setResult("Label", lastRow, "GUV_diameter");
setResult("Mean(membrane)", lastRow, GUV_diameter);
setResult("StdDev(membrane)", lastRow, "%NBD_change");
setResult("Mean(background_new)", lastRow, Ratio);
setOption("ShowRowNumbers", false);
updateResults;

saveAs("Results", MacroDir+OrgFileName+"_Results.csv");
roiManager("Deselect");
roiManager('save', MacroDir+OrgFileName+"_RoiSet.zip");
roiManager("Delete");
run("Close All");
selectWindow("Results");
run("Close");
}

// Generation of log file with data for all files from the analysed folder
selectWindow("Log");
saveAs("Text", dir+"Log.txt");
selectWindow("Log");
run("Close");
