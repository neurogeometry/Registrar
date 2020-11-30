/*
 * MATLAB Compiler: 7.0.1 (R2019a)
 * Date: Thu Sep 19 13:58:50 2019
 * Arguments: 
 * "-B""macro_default""-W""java:StackOverlaps,Class1""-T""link:lib""-d""C:\\Users\\Seyed\\Documents\\DatasetTests\\registrar\\registrar\\JavaMapReduce\\StackOverlaps\\for_testing""class{Class1:C:\\Users\\Seyed\\Documents\\DatasetTests\\registrar\\registrar\\JavaMapReduce\\StackOverlaps.m}"
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
    private static final String sComponentId = "StackOverlap_2E9C11ED5264BCDB600C547915904D4A";
    
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
            new int[]{9,6,0}
        );
    }
    
    public static MWMCR newInstance() throws MWException
    {
        return newInstance(sDefaultComponentOptions);
    }
}
