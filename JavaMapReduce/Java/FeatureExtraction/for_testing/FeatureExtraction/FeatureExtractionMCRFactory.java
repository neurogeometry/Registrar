/*
 * MATLAB Compiler: 6.5 (R2017b)
 * Date: Mon Jan 28 15:08:09 2019
 * Arguments: 
 * "-B""macro_default""-W""java:FeatureExtraction,Class1""-T""link:lib""-d""C:\\Users\\Seyed\\Documents\\DatasetTests\\registrar\\registrar\\JavaMapReduce\\Java\\FeatureExtraction\\for_testing""class{Class1:C:\\Users\\Seyed\\Documents\\DatasetTests\\registrar\\registrar\\JavaMapReduce\\FeatureExtraction.m}"
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
    private static final String sComponentId = "FeatureExtra_E9C7A81A0D1003891123D66E380D8A56";
    
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
            new int[]{9,3,0}
        );
    }
    
    public static MWMCR newInstance() throws MWException
    {
        return newInstance(sDefaultComponentOptions);
    }
}
