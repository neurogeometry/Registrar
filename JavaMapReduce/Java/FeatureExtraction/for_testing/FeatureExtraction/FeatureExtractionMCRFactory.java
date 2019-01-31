/*
 * MATLAB Compiler: 7.0 (R2018b)
 * Date: Thu Jan 31 15:36:44 2019
 * Arguments: 
 * "-B""macro_default""-W""java:FeatureExtraction,Class1""-T""link:lib""-d""/Users/chuhan/Documents/MATLAB/FeatureExtraction/for_testing""class{Class1:/Users/chuhan/Documents/study/2019Spring/registrar/JavaMapReduce/FeatureExtraction.m}"
 */

package FeatureExtraction;

import com.mathworks.toolbox.javabuilder.*;
import com.mathworks.toolbox.javabuilder.internal.*;

/**
 * <i>INTERNAL USE ONLY</i>
 */
public class FeatureExtractionMCRFactory
{
   
    
    /** Component's uuid */
    private static final String sComponentId = "FeatureExtra_69F582BB2D99D7A325102A27F7B9003B";
    
    /** Component name */
    private static final String sComponentName = "FeatureExtraction";
    
   
    /** Pointer to default component options */
    private static final MWComponentOptions sDefaultComponentOptions = 
        new MWComponentOptions(
            MWCtfExtractLocation.EXTRACT_TO_CACHE, 
            new MWCtfClassLoaderSource(FeatureExtractionMCRFactory.class)
        );
    
    
    private FeatureExtractionMCRFactory()
    {
        // Never called.
    }
    
    public static MWMCR newInstance(MWComponentOptions componentOptions) throws MWException
    {
        if (null == componentOptions.getCtfSource()) {
            componentOptions = new MWComponentOptions(componentOptions);
            componentOptions.setCtfSource(sDefaultComponentOptions.getCtfSource());
        }
        return MWMCR.newInstance(
            componentOptions, 
            FeatureExtractionMCRFactory.class, 
            sComponentName, 
            sComponentId,
            new int[]{9,5,0}
        );
    }
    
    public static MWMCR newInstance() throws MWException
    {
        return newInstance(sDefaultComponentOptions);
    }
}
