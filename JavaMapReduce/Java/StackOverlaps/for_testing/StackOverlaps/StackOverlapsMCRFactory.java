/*
 * MATLAB Compiler: 7.0 (R2018b)
 * Date: Thu Jan 31 16:15:54 2019
 * Arguments: 
 * "-B""macro_default""-W""java:StackOverlaps,Class1""-T""link:lib""-d""/Users/chuhan/Documents/MATLAB/StackOverlaps/for_testing""class{Class1:/Users/chuhan/Documents/study/2019Spring/registrar/JavaMapReduce/StackOverlaps.m}"
 */

package StackOverlaps;

import com.mathworks.toolbox.javabuilder.*;
import com.mathworks.toolbox.javabuilder.internal.*;

/**
 * <i>INTERNAL USE ONLY</i>
 */
public class StackOverlapsMCRFactory
{
   
    
    /** Component's uuid */
    private static final String sComponentId = "StackOverlap_4DDCA61451D97EC237519D55E35C7CC4";
    
    /** Component name */
    private static final String sComponentName = "StackOverlaps";
    
   
    /** Pointer to default component options */
    private static final MWComponentOptions sDefaultComponentOptions = 
        new MWComponentOptions(
            MWCtfExtractLocation.EXTRACT_TO_CACHE, 
            new MWCtfClassLoaderSource(StackOverlapsMCRFactory.class)
        );
    
    
    private StackOverlapsMCRFactory()
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
            StackOverlapsMCRFactory.class, 
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
