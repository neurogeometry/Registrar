/*
 * MATLAB Compiler: 7.0 (R2018b)
 * Date: Wed Feb 20 18:26:27 2019
 * Arguments: 
 * "-B""macro_default""-W""java:FeatureExtractionFunc,Class3""-T""link:lib""-d""/Users/chuhan/Documents/study/2019Spring/registrar/JavaMapReduce/FeatureExtractionFunc/for_testing""class{Class3:/Users/chuhan/Documents/study/2019Spring/registrar/JavaMapReduce/FeatureExtractionFunc.m}"
 */

package FeatureExtractionFunc;

import com.mathworks.toolbox.javabuilder.*;
import com.mathworks.toolbox.javabuilder.internal.*;

/**
 * <i>INTERNAL USE ONLY</i>
 */
public class FeatureExtractionFuncMCRFactory
{
   
    
    /** Component's uuid */
    private static final String sComponentId = "FeatureExtra_6993DB45EDB82368F40238D1ABEC83B1";
    
    /** Component name */
    private static final String sComponentName = "FeatureExtractionFunc";
    
   
    /** Pointer to default component options */
    private static final MWComponentOptions sDefaultComponentOptions = 
        new MWComponentOptions(
            MWCtfExtractLocation.EXTRACT_TO_CACHE, 
            new MWCtfClassLoaderSource(FeatureExtractionFuncMCRFactory.class)
        );
    
    
    private FeatureExtractionFuncMCRFactory()
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
            FeatureExtractionFuncMCRFactory.class, 
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
