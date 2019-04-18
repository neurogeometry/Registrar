/*
 * MATLAB Compiler: 7.0 (R2018b)
 * Date: Thu Apr 18 14:39:35 2019
 * Arguments: 
 * "-B""macro_default""-W""java:DataRetiling,MATDataRetiling""-T""link:lib""-d""/Users/anubh/Documents/NCTracer/registrar/JavaMapReduce/DataRetiling/for_testing""class{MATDataRetiling:/Users/anubh/Documents/NCTracer/registrar/JavaMapReduce/DataRetiling.m}"
 */

package DataRetiling;

import com.mathworks.toolbox.javabuilder.*;
import com.mathworks.toolbox.javabuilder.internal.*;

/**
 * <i>INTERNAL USE ONLY</i>
 */
public class DataRetilingMCRFactory
{
   
    
    /** Component's uuid */
    private static final String sComponentId = "DataRetiling_1DC5244D395DDAC490BF914B33870B55";
    
    /** Component name */
    private static final String sComponentName = "DataRetiling";
    
   
    /** Pointer to default component options */
    private static final MWComponentOptions sDefaultComponentOptions = 
        new MWComponentOptions(
            MWCtfExtractLocation.EXTRACT_TO_CACHE, 
            new MWCtfClassLoaderSource(DataRetilingMCRFactory.class)
        );
    
    
    private DataRetilingMCRFactory()
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
            DataRetilingMCRFactory.class, 
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
