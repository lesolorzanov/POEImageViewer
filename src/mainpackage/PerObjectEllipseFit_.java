/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package mainpackage;

import java.awt.*;
import java.awt.color.*;
import java.awt.image.*;
import java.io.*;
import javax.imageio.ImageIO;

/**
 * This plugin segments objects using local adaptive thresholding optimizing
 * Ellipsefit Ranefall P, Sadanandan SK, Wählby C. "Fast Adaptive Local
 * Thresholding Based on Ellipse Fit", (International Symposium on Biomedical
 * Imaging (ISBI'16), Prague, Czech Republic, April 13-16, 2016)
 *
 * @author Petter Ranefall - Leslie Solorzano
 * @version 1.0
 * @date 2017-4-27
 *
 * Input is an interval of expected object sizes. The algorithm computes the per
 * object threshold level that gives the maximum ellipse fit. Additional
 * parameters are allowed intervals of major- and minor axes, the ratio between
 * those, and lowest ellipse fit.
 */
public class PerObjectEllipseFit_ {

    /*protected ImageStack stack;
    protected ImagePlus imp;*/
    
    BufferedImage image;
    BufferedImage maskimage;
    
    byte state;
    /*int width;
    int height;*/

    int[] parNode;
    short[] imArray;
    byte[] outputArray;
    double[] ellipseFit;
    double[] majorAxis;
    double[] minorAxis;
    short maxP = 0, minP = 0;
    int minSize = 0;
    int maxSize = 0;
    int minMajorAxis = 0;
    int maxMajorAxis = 0;
    int minMinorAxis = 0;
    int maxMinorAxis = 0;
    double minMajorMinorRatio = 1.0;
    double maxMajorMinorRatio = 10.0;
    
    int nPixels;

    double ellipseThr = 0;
    int[] areas;
    int[] sortedIndex;
    double[] sumX;
    double[] sumX2;
    double[] sumY;
    double[] sumY2;
    double[] sumXY;
    
    int width;
    int height;
    int nSlices;
    boolean darkBkg = false;

    //CONSTANTS
    public static final byte UNDEFINED = (byte) 0;
    public static final byte NOT_OBJECT = (byte) 1;
    public static final byte MAYBE_NOT_OBJECT = (byte) 2;
    public static final byte OBJECT = (byte) 3;

    public PerObjectEllipseFit_(File f) {
        
        try{
            image = ImageIO.read(f);
            width = image.getWidth();
            height = image.getHeight(); 
            nPixels = width*height;
            state=1; //image is loaded
        }catch(Exception e){
            System.out.println(e.getStackTrace());
        }
        
    }
    
    public PerObjectEllipseFit_() {
        
        state=0; //image is not loaded
        
    }
    
    public void init(String name){
        try{
            System.out.println(name);
            setImage(ImageIO.read(new File(name)));
            width = image.getWidth();
            height = image.getHeight(); 
            nPixels = width*height;
        }catch(Exception e){
            System.out.println(e.getStackTrace());
        }
        
        state=1; //image is loaded
        createTree();
        
        
        
    }

    
    public void updateVars(int minSize, int maxSize, double ellipseThr, 
        int minMajorAxis, int maxMajorAxis, int minMinorAxis, int maxMinorAxis, 
        double minMajorMinorRatio, double maxMajorMinorRatio, boolean darkBkg ){
        
        setMinSize(minSize);
        setMaxSize(maxSize);
        setEllipseThr(ellipseThr);
        setMinMajorAxis(minMajorAxis);
        setMaxMajorAxis(maxMajorAxis);
        setMinMinorAxis(minMinorAxis);
        setMaxMinorAxis(maxMinorAxis);
        setMinMajorMinorRatio(minMajorMinorRatio);
        setMaxMajorMinorRatio(maxMajorMinorRatio);
        setDarkBkg(darkBkg);
    
    }
    
 
    public void createTree(){
        
        //System.out.print("fbgdfgdfgdf"); 
         //Create nodes
        parNode = new int[nPixels];

        // 1) -----------------------------
        //int type = imp.getType();
        //int slice = 1;
        int x, y, x0, x2, y0, y2;
        //int rest;
        int i, j, k;
        short p;
        int q;
        imArray = new short[nPixels];
        outputArray = new byte[nPixels];
        ellipseFit = new double[nPixels];
        majorAxis = new double[nPixels];
        minorAxis = new double[nPixels];
        sumX = new double[nPixels];
        sumX2 = new double[nPixels];
        sumY = new double[nPixels];
        sumY2 = new double[nPixels];
        sumXY = new double[nPixels];
        
        // 2) -----------------------------
        
        int mask = 0xff;
        //Find min, max. Copy to imArray
        //byte[] pixels = (byte[]) stack.getPixels(1);
        byte[] pixels = ((DataBufferByte)image.getRaster().getDataBuffer()).getData() ;
        p = (short) (mask & pixels[0]);
        maxP = minP = p;
        //pixels = (byte[]) stack.getPixels(1);
        for (i = 0; i < nPixels; i++) {
            p = (short) (mask & pixels[i]);
            if (p < minP) {
                minP = p;
            }
            if (p > maxP) {
                maxP = p;
            }
            imArray[i] = p;
        }
        
        /*for (int l = 0; l < 20; l++) {
           System.out.print(imArray[l]+" ");
        }*/
        
        // 3) -----------------------------
        if (!darkBkg) {
            for (i = 0; i < nPixels; i++) {
                imArray[i] = (short) ((int) minP + (int) maxP - (int) imArray[i]);
            }
        }
        // 4) -----------------------------
        areas = new int[nPixels];
        for (i = 0; i < nPixels; i++) {
            areas[i] = 1;
        }
        i = 0;
        // 5) -----------------------------
        for (y = 0; y < height; y++) {
            for (x = 0; x < width; x++) {
                sumX[i] = x;
                sumX2[i] = x * x;
                sumY[i] = y;
                sumY2[i] = y * y;
                sumXY[i] = x * y;
                i++;
            }
        }
        // 6) -----------------------------
        //Sort points
        // create a counting array, counts, with a member for
        // each possible discrete value in the input.
        // initialize all counts to 0.
        int nLevels = maxP - minP + 1;
        int[] counts = new int[nLevels];
        // 6A) -----------------------------
        // for each value in the unsorted array, increment the
        // count in the corresponding element of the count array
        for (i = 0; i < nPixels; i++) {
            counts[imArray[i] - minP]++;
        }
        // 6B) -----------------------------
        // accumulate the counts - the result is that counts will hold
        // the offset into the sorted array for the value associated with that index
        for (i = 1; i < nLevels; i++) {
            counts[i] += counts[i - 1];
        }
        // store the elements in a new ordered array
        sortedIndex = new int[nPixels];
        // 6C) -----------------------------
        for (i = nPixels - 1; i >= 0; i--) {
            // decrementing the counts value ensures duplicate values in A
            // are stored at different indices in sorted.
            sortedIndex[--counts[imArray[i] - minP]] = i;
        }

        // 7) -----------------------------
        //Init nodes
        for (i = 0; i < nPixels; i++) {
            parNode[i] = i;
        }
        //Search in decreasing order
        int curNode;
        int adjNode;
        for (i = nPixels - 1; i >= 0; i--) {
            j = sortedIndex[i];
            curNode = j;

            y = j / width;
            x = j - y * width;

            y0 = y - 1;
            y2 = y + 1;
            x0 = x - 1;
            x2 = x + 1;

            //Later neigbours x2,y2
            if (y2 < height) {
                k = x + width * y2;
                if (imArray[k] >= imArray[j]) {
                    adjNode = findNode(k);
                    if (curNode != adjNode) {
                        curNode = mergeNodes(adjNode, curNode);
                    }
                }
            }
            if (x2 < width) {
                k = x2 + width * y;
                if (imArray[k] >= imArray[j]) {
                    adjNode = findNode(k);
                    if (curNode != adjNode) {
                        curNode = mergeNodes(adjNode, curNode);
                    }
                }
            }
            //Earlier neighbours x0,y0. No need to check =
            if (x0 >= 0) {
                k = x0 + width * y;
                if (imArray[k] > imArray[j]) {
                    adjNode = findNode(k);
                    if (curNode != adjNode) {
                        curNode = mergeNodes(adjNode, curNode);
                    }

                }
            }
            if (y0 >= 0) {
                k = x + width * y0;
                if (imArray[k] > imArray[j]) {
                    adjNode = findNode(k);
                    if (curNode != adjNode) {
                        curNode = mergeNodes(adjNode, curNode);
                    }
                }
            }

        }
        
        state=2; //tree is created
    
    }   
    
    public void exec() {
        
        if(state<2){
            System.out.println("Tree is not created");
            return;            
            //if tree is not created we cannot execute            
        } //tree is created
        
        //int rest;
        int i, j;
       
        for (i = nPixels - 1; i >= 0; i--) {
            j = sortedIndex[i];
            if (outputArray[j] == UNDEFINED) {
                int e = j;
                while (imArray[e] == imArray[parNode[e]] && outputArray[e] == UNDEFINED) {
                    e = parNode[e];
                }
                if (outputArray[e] == UNDEFINED) {
                    findBestEllipseLevel(j, minSize, maxSize, ellipseThr, minMajorAxis, maxMajorAxis, minMinorAxis, maxMinorAxis, minMajorMinorRatio, maxMajorMinorRatio, -1);
                } else {
                    int e1 = j;
                    while (e1 != e) {
                        outputArray[e1] = outputArray[e];
                        e1 = parNode[e1];
                    }
                }
            }
        }
        //Handle MAYBE
        for (i = nPixels - 1; i >= 0; i--) {
            j = sortedIndex[i];
            if (outputArray[j] == MAYBE_NOT_OBJECT) {
                findBestEllipseLevel(j, minSize, maxSize, ellipseThr, minMajorAxis, maxMajorAxis, minMinorAxis, maxMinorAxis, minMajorMinorRatio, maxMajorMinorRatio, -1);
            }
        }
        for (i = nPixels - 1; i >= 0; i--) {
            j = sortedIndex[i];
            if (outputArray[j] == MAYBE_NOT_OBJECT) {
                findObjBelow(j, maxSize);
            }
        }
        
        for (i=0; i<nPixels; i++)
        {
            if(outputArray[i] == OBJECT)
                outputArray[i]=(byte)255;
            else
                outputArray[i]=(byte)0;
        }

        //create the mask image
        
        //maskimage = new BufferedImage(width, height,BufferedImage.TYPE_BYTE_BINARY);
        
        
        ColorSpace cs = ColorSpace.getInstance(ColorSpace.CS_GRAY);
        int[] nBits = { 8 };
        ColorModel cm = new ComponentColorModel(cs, nBits, false, true, Transparency.OPAQUE, DataBuffer.TYPE_BYTE);
        SampleModel sm = cm.createCompatibleSampleModel(width, height);
        DataBufferByte db = new DataBufferByte(outputArray, width * height);
        WritableRaster raster = Raster.createWritableRaster(sm, db, null);
        maskimage = new BufferedImage(cm, raster, false, null);

        
        /*final byte[] a = ( (DataBufferByte) maskimage.getRaster().getDataBuffer() ).getData();
        System.arraycopy(outputArray, 0, a, 0, outputArray.length);*/
        
        /*for (int x = 0; x < 10; x++) {
            for (int y = 0; y < 10; y++) {
                maskimage.setRGB(x, y, outputArray[x+y*width]);
            }            
        }*/
           
    }
    
    public void setOutputArrayToZero(){
        if(state>=2){
            for (int i = 0; i < outputArray.length; i++) {
                outputArray[i]=UNDEFINED;
            }
        }
    
    }

    public int findNode(int e) {
        if (parNode[e] != e) {
            int root = findNode(parNode[e]);
            //parNode[e] = root; //This cannot be used here
            return root;
        } else {
            return e;
        }
    }

    public boolean findObjBelow(int e, double maxS) {
        int y = e / width;
        int x = e - y * width;
        if (parNode[e] == e) {
            outputArray[e] = NOT_OBJECT;
            return false;
        }
        if (outputArray[parNode[e]] == OBJECT) {
            outputArray[e] = OBJECT;
            return true;
        }
        if (areas[parNode[e]] > maxS) {
            outputArray[e] = NOT_OBJECT;
            return false;
        }
        boolean found = findObjBelow(parNode[e], maxS);
        if (found) {
            outputArray[e] = outputArray[parNode[e]];
        } else {
            outputArray[e] = NOT_OBJECT;
        }
        return found;
    }

    public int findBestEllipseLevel(int e, double minS, double maxS, double eThr, 
            double minMajor, double maxMajor, double minMinor, 
            double maxMinor, double minMajorMinor, double maxMajorMinor, int optE) {
        double optEf;
        double optArea;
        if (optE >= 0) {
            optEf = ellipseFit[optE];
            optArea = areas[optE];
        } else {
            optEf = 0;
            optArea = 0;
        }
        int startE = e;
        while (imArray[e] == imArray[parNode[e]] && parNode[e] != e) {
            e = parNode[e];
        }
        if (outputArray[e] == OBJECT) {
            return e;
        }
        if (outputArray[e] == NOT_OBJECT) {
            return optE;
        }
        if (parNode[e] == e) {
            outputArray[e] = NOT_OBJECT;
            return optE;
        }
        if (areas[e] > maxS) {
            outputArray[startE] = NOT_OBJECT;
            while (startE != parNode[startE] && outputArray[parNode[startE]] == UNDEFINED) {
                startE = parNode[startE];
                outputArray[startE] = NOT_OBJECT;
            }
            return optE;
        }
        double majorMinorRatio = majorAxis[e] / minorAxis[e];
        if (areas[e] >= minS && ellipseFit[e] > optEf && ellipseFit[e] > eThr && majorAxis[e] >= minMajor && majorAxis[e] <= maxMajor && minorAxis[e] >= minMinor && minorAxis[e] <= maxMinor && majorMinorRatio >= minMajorMinor && majorMinorRatio <= maxMajorMinor) {
            optEf = ellipseFit[e];
            optArea = areas[e];
            optE = e;
        }
        optE = findBestEllipseLevel(parNode[e], minS, maxS, eThr, minMajor, maxMajor, minMinor, maxMinor, minMajorMinor, maxMajorMinor, optE);
        if (optE >= 0) {
            optEf = ellipseFit[optE];
            optArea = areas[optE];
        } else {
            optEf = 0;
            optArea = 0;
        }
        if (optE >= 0 && imArray[e] >= imArray[optE]) //Opt found below
        {
            while (startE != e) {
                outputArray[startE] = OBJECT;
                startE = parNode[startE];
            }
            outputArray[e] = OBJECT;
            return optE;
        }

        // Below best level
        outputArray[startE] = MAYBE_NOT_OBJECT;
        if (optArea > areas[startE]) {
            areas[startE] = 0;
        } else {
            areas[startE] -= optArea;
        }

        while (startE != e) {
            startE = parNode[startE];
            outputArray[startE] = MAYBE_NOT_OBJECT;
            if (optArea > areas[startE]) {
                areas[startE] = 0;
            } else {
                areas[startE] -= optArea;
            }
        }
        return optE;
    }

    public int mergeNodes(int e1, int e2) {
        int res;
        int m;

        if (imArray[e1] == imArray[e2]) {
            res = Math.max(e1, e2);
            m = Math.min(e1, e2);
        } else {
            res = e2;
            m = e1;
        }
        areas[res] += areas[m];
        parNode[m] = res;
        sumX[res] += sumX[m];
        sumX2[res] += sumX2[m];
        sumY[res] += sumY[m];
        sumY2[res] += sumY2[m];
        sumXY[res] += sumXY[m];
        if (areas[res] > 1) {
            double varX = (sumX2[res] - sumX[res] * sumX[res] / areas[res]);
            double varY = (sumY2[res] - sumY[res] * sumY[res] / areas[res]);
            double covXY = (sumXY[res] - sumX[res] * sumY[res] / areas[res]);
            double varXPlusVarY = varX + varY;
            double sqrtExpr = Math.sqrt(varXPlusVarY * varXPlusVarY - 4 * (varX * varY - covXY * covXY));
            if (varXPlusVarY > sqrtExpr) {
                double r1 = 2.0 * Math.sqrt((varXPlusVarY + sqrtExpr) / (2 * areas[res]));
                double r2 = 2.0 * Math.sqrt((varXPlusVarY - sqrtExpr) / (2 * areas[res]));
                double ellipseArea = Math.PI * r1 * r2;
                if (ellipseArea > 0) {
                    ellipseFit[res] = areas[res] / ellipseArea;
                    majorAxis[res] = 2.0 * r1;
                    minorAxis[res] = 2.0 * r2;
                }
            }
        }
        return res;
    }

    public BufferedImage getImage() {
        return image;
    }

    public void setImage(BufferedImage image) {
        this.image = image;
    }

    public BufferedImage getMaskimage() {
        return maskimage;
    }

    public void setMaskimage(BufferedImage maskimage) {
        this.maskimage = maskimage;
    }
     
    public int getMinSize() {
        return minSize;
    }

    public void setMinSize(int minSize) {
        this.minSize = minSize;
    }

    public int getMaxSize() {
        return maxSize;
    }

    public void setMaxSize(int maxSize) {
        this.maxSize = maxSize;
    }

    public int getMinMajorAxis() {
        return minMajorAxis;
    }

    public void setMinMajorAxis(int minMajorAxis) {
        this.minMajorAxis = minMajorAxis;
    }

    public int getMaxMajorAxis() {
        return maxMajorAxis;
    }

    public void setMaxMajorAxis(int maxMajorAxis) {
        this.maxMajorAxis = maxMajorAxis;
    }

    public int getMinMinorAxis() {
        return minMinorAxis;
    }

    public void setMinMinorAxis(int minMinorAxis) {
        this.minMinorAxis = minMinorAxis;
    }

    public int getMaxMinorAxis() {
        return maxMinorAxis;
    }

    public void setMaxMinorAxis(int maxMinorAxis) {
        this.maxMinorAxis = maxMinorAxis;
    }

    public double getMinMajorMinorRatio() {
        return minMajorMinorRatio;
    }

    public void setMinMajorMinorRatio(double minMajorMinorRatio) {
        this.minMajorMinorRatio = minMajorMinorRatio;
    }

    public double getMaxMajorMinorRatio() {
        return maxMajorMinorRatio;
    }

    public void setMaxMajorMinorRatio(double maxMajorMinorRatio) {
        this.maxMajorMinorRatio = maxMajorMinorRatio;
    }

    public double getEllipseThr() {
        return ellipseThr;
    }

    public void setEllipseThr(double ellipseThr) {
        this.ellipseThr = ellipseThr;
    }

    public boolean isDarkBkg() {
        return darkBkg;
    }

    public void setDarkBkg(boolean darkBkg) {
        this.darkBkg = darkBkg;
    }
    
    
    
}
