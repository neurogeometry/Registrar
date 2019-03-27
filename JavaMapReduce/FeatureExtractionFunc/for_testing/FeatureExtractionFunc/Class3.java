/*
 * MATLAB Compiler: 7.0 (R2018b)
 * Date: Wed Feb 20 18:26:27 2019
 * Arguments: 
 * "-B""macro_default""-W""java:FeatureExtractionFunc,Class3""-T""link:lib""-d""/Users/chuhan/Documents/study/2019Spring/registrar/JavaMapReduce/FeatureExtractionFunc/for_testing""class{Class3:/Users/chuhan/Documents/study/2019Spring/registrar/JavaMapReduce/FeatureExtractionFunc.m}"
 */

package FeatureExtractionFunc;

import com.mathworks.toolbox.javabuilder.*;
import com.mathworks.toolbox.javabuilder.internal.*;
import java.util.*;

/**
 * The <code>Class3</code> class provides a Java interface to MATLAB functions. 
 * The interface is compiled from the following files:
 * <pre>
 *  /Users/chuhan/Documents/study/2019Spring/registrar/JavaMapReduce/FeatureExtractionFunc.m
 * </pre>
 * The {@link #dispose} method <b>must</b> be called on a <code>Class3</code> instance 
 * when it is no longer needed to ensure that native resources allocated by this class 
 * are properly freed.
 * @version 0.0
 */
public class Class3 extends MWComponentInstance<Class3>
{
    /**
     * Tracks all instances of this class to ensure their dispose method is
     * called on shutdown.
     */
    private static final Set<Disposable> sInstances = new HashSet<Disposable>();

    /**
     * Maintains information used in calling the <code>FeatureExtractionFunc</code> 
     *MATLAB function.
     */
    private static final MWFunctionSignature sFeatureExtractionFuncSignature =
        new MWFunctionSignature(/* max outputs = */ 1,
                                /* has varargout = */ false,
                                /* function name = */ "FeatureExtractionFunc",
                                /* max inputs = */ 5,
                                /* has varargin = */ false);

    /**
     * Shared initialization implementation - private
     * @throws MWException An error has occurred during the function call.
     */
    private Class3 (final MWMCR mcr) throws MWException
    {
        super(mcr);
        // add this to sInstances
        synchronized(Class3.class) {
            sInstances.add(this);
        }
    }

    /**
     * Constructs a new instance of the <code>Class3</code> class.
     * @throws MWException An error has occurred during the function call.
     */
    public Class3() throws MWException
    {
        this(FeatureExtractionFuncMCRFactory.newInstance());
    }
    
    private static MWComponentOptions getPathToComponentOptions(String path)
    {
        MWComponentOptions options = new MWComponentOptions(new MWCtfExtractLocation(path),
                                                            new MWCtfDirectorySource(path));
        return options;
    }
    
    /**
     * @deprecated Please use the constructor {@link #Class3(MWComponentOptions componentOptions)}.
     * The <code>com.mathworks.toolbox.javabuilder.MWComponentOptions</code> class provides an API to set the
     * path to the component.
     * @param pathToComponent Path to component directory.
     * @throws MWException An error has occurred during the function call.
     */
    public Class3(String pathToComponent) throws MWException
    {
        this(FeatureExtractionFuncMCRFactory.newInstance(getPathToComponentOptions(pathToComponent)));
    }
    
    /**
     * Constructs a new instance of the <code>Class3</code> class. Use this constructor 
     * to specify the options required to instantiate this component.  The options will 
     * be specific to the instance of this component being created.
     * @param componentOptions Options specific to the component.
     * @throws MWException An error has occurred during the function call.
     */
    public Class3(MWComponentOptions componentOptions) throws MWException
    {
        this(FeatureExtractionFuncMCRFactory.newInstance(componentOptions));
    }
    
    /** Frees native resources associated with this object */
    public void dispose()
    {
        try {
            super.dispose();
        } finally {
            synchronized(Class3.class) {
                sInstances.remove(this);
            }
        }
    }
  
    /**
     * Invokes the first MATLAB function specified to MCC, with any arguments given on
     * the command line, and prints the result.
     *
     * @param args arguments to the function
     */
    public static void main (String[] args)
    {
        try {
            MWMCR mcr = FeatureExtractionFuncMCRFactory.newInstance();
            mcr.runMain( sFeatureExtractionFuncSignature, args);
            mcr.dispose();
        } catch (Throwable t) {
            t.printStackTrace();
        }
    }
    
    /**
     * Calls dispose method for each outstanding instance of this class.
     */
    public static void disposeAllInstances()
    {
        synchronized(Class3.class) {
            for (Disposable i : sInstances) i.dispose();
            sInstances.clear();
        }
    }

    /**
     * Provides the interface for calling the <code>FeatureExtractionFunc</code> MATLAB function 
     * where the first argument, an instance of List, receives the output of the MATLAB function and
     * the second argument, also an instance of List, provides the input to the MATLAB function.
     * <p>
     * Description as provided by the author of the MATLAB function:
     * </p>
     * <pre>
     * % -------------------------------------------------------------------------
     * % Author: Seyed Mostafa Mousavi Kahaki, Armen Stepanyants
     * % Northeastern University, USA
     * % =========================================================================
     * % -------------------------------------------------------------------------
     * </pre>
     * @param lhs List in which to return outputs. Number of outputs (nargout) is
     * determined by allocated size of this List. Outputs are returned as
     * sub-classes of <code>com.mathworks.toolbox.javabuilder.MWArray</code>.
     * Each output array should be freed by calling its <code>dispose()</code>
     * method.
     *
     * @param rhs List containing inputs. Number of inputs (nargin) is determined
     * by the allocated size of this List. Input arguments may be passed as
     * sub-classes of <code>com.mathworks.toolbox.javabuilder.MWArray</code>, or
     * as arrays of any supported Java type. Arguments passed as Java types are
     * converted to MATLAB arrays according to default conversion rules.
     * @throws MWException An error has occurred during the function call.
     */
    public void FeatureExtractionFunc(List lhs, List rhs) throws MWException
    {
        fMCR.invoke(lhs, rhs, sFeatureExtractionFuncSignature);
    }

    /**
     * Provides the interface for calling the <code>FeatureExtractionFunc</code> MATLAB function 
     * where the first argument, an Object array, receives the output of the MATLAB function and
     * the second argument, also an Object array, provides the input to the MATLAB function.
     * <p>
     * Description as provided by the author of the MATLAB function:
     * </p>
     * <pre>
     * % -------------------------------------------------------------------------
     * % Author: Seyed Mostafa Mousavi Kahaki, Armen Stepanyants
     * % Northeastern University, USA
     * % =========================================================================
     * % -------------------------------------------------------------------------
     * </pre>
     * @param lhs array in which to return outputs. Number of outputs (nargout)
     * is determined by allocated size of this array. Outputs are returned as
     * sub-classes of <code>com.mathworks.toolbox.javabuilder.MWArray</code>.
     * Each output array should be freed by calling its <code>dispose()</code>
     * method.
     *
     * @param rhs array containing inputs. Number of inputs (nargin) is
     * determined by the allocated size of this array. Input arguments may be
     * passed as sub-classes of
     * <code>com.mathworks.toolbox.javabuilder.MWArray</code>, or as arrays of
     * any supported Java type. Arguments passed as Java types are converted to
     * MATLAB arrays according to default conversion rules.
     * @throws MWException An error has occurred during the function call.
     */
    public void FeatureExtractionFunc(Object[] lhs, Object[] rhs) throws MWException
    {
        fMCR.invoke(Arrays.asList(lhs), Arrays.asList(rhs), sFeatureExtractionFuncSignature);
    }

    /**
     * Provides the standard interface for calling the <code>FeatureExtractionFunc</code> MATLAB function with 
     * 5 comma-separated input arguments.
     * Input arguments may be passed as sub-classes of
     * <code>com.mathworks.toolbox.javabuilder.MWArray</code>, or as arrays of
     * any supported Java type. Arguments passed as Java types are converted to
     * MATLAB arrays according to default conversion rules.
     *
     * <p>
     * Description as provided by the author of the MATLAB function:
     * </p>
     * <pre>
     * % -------------------------------------------------------------------------
     * % Author: Seyed Mostafa Mousavi Kahaki, Armen Stepanyants
     * % Northeastern University, USA
     * % =========================================================================
     * % -------------------------------------------------------------------------
     * </pre>
     * @param nargout Number of outputs to return.
     * @param rhs The inputs to the MATLAB function.
     * @return Array of length nargout containing the function outputs. Outputs
     * are returned as sub-classes of
     * <code>com.mathworks.toolbox.javabuilder.MWArray</code>. Each output array
     * should be freed by calling its <code>dispose()</code> method.
     * @throws MWException An error has occurred during the function call.
     */
    public Object[] FeatureExtractionFunc(int nargout, Object... rhs) throws MWException
    {
        Object[] lhs = new Object[nargout];
        fMCR.invoke(Arrays.asList(lhs), 
                    MWMCR.getRhsCompat(rhs, sFeatureExtractionFuncSignature), 
                    sFeatureExtractionFuncSignature);
        return lhs;
    }
}
