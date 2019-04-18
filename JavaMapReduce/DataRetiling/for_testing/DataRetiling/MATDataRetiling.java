/*
 * MATLAB Compiler: 7.0 (R2018b)
 * Date: Thu Apr 18 14:39:35 2019
 * Arguments: 
 * "-B""macro_default""-W""java:DataRetiling,MATDataRetiling""-T""link:lib""-d""/Users/anubh/Documents/NCTracer/registrar/JavaMapReduce/DataRetiling/for_testing""class{MATDataRetiling:/Users/anubh/Documents/NCTracer/registrar/JavaMapReduce/DataRetiling.m}"
 */

package DataRetiling;

import com.mathworks.toolbox.javabuilder.*;
import com.mathworks.toolbox.javabuilder.internal.*;
import java.util.*;

/**
 * The <code>MATDataRetiling</code> class provides a Java interface to MATLAB functions. 
 * The interface is compiled from the following files:
 * <pre>
 *  /Users/anubh/Documents/NCTracer/registrar/JavaMapReduce/DataRetiling.m
 * </pre>
 * The {@link #dispose} method <b>must</b> be called on a <code>MATDataRetiling</code> 
 * instance when it is no longer needed to ensure that native resources allocated by this 
 * class are properly freed.
 * @version 0.0
 */
public class MATDataRetiling extends MWComponentInstance<MATDataRetiling>
{
    /**
     * Tracks all instances of this class to ensure their dispose method is
     * called on shutdown.
     */
    private static final Set<Disposable> sInstances = new HashSet<Disposable>();

    /**
     * Maintains information used in calling the <code>DataRetiling</code> MATLAB 
     *function.
     */
    private static final MWFunctionSignature sDataRetilingSignature =
        new MWFunctionSignature(/* max outputs = */ 0,
                                /* has varargout = */ false,
                                /* function name = */ "DataRetiling",
                                /* max inputs = */ 13,
                                /* has varargin = */ false);

    /**
     * Shared initialization implementation - private
     * @throws MWException An error has occurred during the function call.
     */
    private MATDataRetiling (final MWMCR mcr) throws MWException
    {
        super(mcr);
        // add this to sInstances
        synchronized(MATDataRetiling.class) {
            sInstances.add(this);
        }
    }

    /**
     * Constructs a new instance of the <code>MATDataRetiling</code> class.
     * @throws MWException An error has occurred during the function call.
     */
    public MATDataRetiling() throws MWException
    {
        this(DataRetilingMCRFactory.newInstance());
    }
    
    private static MWComponentOptions getPathToComponentOptions(String path)
    {
        MWComponentOptions options = new MWComponentOptions(new MWCtfExtractLocation(path),
                                                            new MWCtfDirectorySource(path));
        return options;
    }
    
    /**
     * @deprecated Please use the constructor {@link #MATDataRetiling(MWComponentOptions componentOptions)}.
     * The <code>com.mathworks.toolbox.javabuilder.MWComponentOptions</code> class provides an API to set the
     * path to the component.
     * @param pathToComponent Path to component directory.
     * @throws MWException An error has occurred during the function call.
     */
    public MATDataRetiling(String pathToComponent) throws MWException
    {
        this(DataRetilingMCRFactory.newInstance(getPathToComponentOptions(pathToComponent)));
    }
    
    /**
     * Constructs a new instance of the <code>MATDataRetiling</code> class. Use this 
     * constructor to specify the options required to instantiate this component.  The 
     * options will be specific to the instance of this component being created.
     * @param componentOptions Options specific to the component.
     * @throws MWException An error has occurred during the function call.
     */
    public MATDataRetiling(MWComponentOptions componentOptions) throws MWException
    {
        this(DataRetilingMCRFactory.newInstance(componentOptions));
    }
    
    /** Frees native resources associated with this object */
    public void dispose()
    {
        try {
            super.dispose();
        } finally {
            synchronized(MATDataRetiling.class) {
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
            MWMCR mcr = DataRetilingMCRFactory.newInstance();
            mcr.runMain( sDataRetilingSignature, args);
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
        synchronized(MATDataRetiling.class) {
            for (Disposable i : sInstances) i.dispose();
            sInstances.clear();
        }
    }

    /**
     * Provides the interface for calling the <code>DataRetiling</code> MATLAB function 
     * where the first argument, an instance of List, receives the output of the MATLAB function and
     * the second argument, also an instance of List, provides the input to the MATLAB function.
     * <p>
     * Description as provided by the author of the MATLAB function:
     * </p>
     * <pre>
     * % paramsRERemoveBlack = 1;
     * % paramsBigTileSize = [512 512 128]; % powers of 2 only
     * % paramsFinalTileSize = [128*4 128*4 128]; % powers of 2 only
     * % paramsEmptyVoxelsValue=111;
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
    public void DataRetiling(List lhs, List rhs) throws MWException
    {
        fMCR.invoke(lhs, rhs, sDataRetilingSignature);
    }

    /**
     * Provides the interface for calling the <code>DataRetiling</code> MATLAB function 
     * where the first argument, an Object array, receives the output of the MATLAB function and
     * the second argument, also an Object array, provides the input to the MATLAB function.
     * <p>
     * Description as provided by the author of the MATLAB function:
     * </p>
     * <pre>
     * % paramsRERemoveBlack = 1;
     * % paramsBigTileSize = [512 512 128]; % powers of 2 only
     * % paramsFinalTileSize = [128*4 128*4 128]; % powers of 2 only
     * % paramsEmptyVoxelsValue=111;
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
    public void DataRetiling(Object[] lhs, Object[] rhs) throws MWException
    {
        fMCR.invoke(Arrays.asList(lhs), Arrays.asList(rhs), sDataRetilingSignature);
    }

    /**
     * Provides the standard interface for calling the <code>DataRetiling</code> MATLAB function with 
     * 13 comma-separated input arguments.
     * Input arguments may be passed as sub-classes of
     * <code>com.mathworks.toolbox.javabuilder.MWArray</code>, or as arrays of
     * any supported Java type. Arguments passed as Java types are converted to
     * MATLAB arrays according to default conversion rules.
     *
     * <p>
     * Description as provided by the author of the MATLAB function:
     * </p>
     * <pre>
     * % paramsRERemoveBlack = 1;
     * % paramsBigTileSize = [512 512 128]; % powers of 2 only
     * % paramsFinalTileSize = [128*4 128*4 128]; % powers of 2 only
     * % paramsEmptyVoxelsValue=111;
     * </pre>
     * @param rhs The inputs to the MATLAB function.
     * @return Array of length nargout containing the function outputs. Outputs
     * are returned as sub-classes of
     * <code>com.mathworks.toolbox.javabuilder.MWArray</code>. Each output array
     * should be freed by calling its <code>dispose()</code> method.
     * @throws MWException An error has occurred during the function call.
     */
    public Object[] DataRetiling(Object... rhs) throws MWException
    {
        Object[] lhs = new Object[0];
        fMCR.invoke(Arrays.asList(lhs), 
                    MWMCR.getRhsCompat(rhs, sDataRetilingSignature), 
                    sDataRetilingSignature);
        return lhs;
    }
}
