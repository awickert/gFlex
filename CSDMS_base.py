
## Copyright (c) 2001-2010, Scott D. Peckham
## August 2009
## May 2010 (initialize_config_vars(), cleanup, etc.)

#-----------------------------------------------------------------------
#  Notes:  This file defines a "base class" for CSDMS "process"
#          components.  Many of these functions are implementations
#          of methods defined in "topoflow3.IRFPort.sidl".  They are 
#          therefore required to use TF components in a CCA framework.

#  Idea:   Rename "get_cca_ports" to "get_cca_uses_ports".
#-----------------------------------------------------------------------
#
#  unit_test()
#  sync_file_test()   ## (5/7/10)  OBSOLETE ???
#
#  class CSDMS_component
#
#      __init__()
#      get_status()
#      is_scalar()
#      is_vector()
#      is_grid()
#      has_variable()
#      -----------------------------
#       Next 3 currently identical
#      -----------------------------
#      get_scalar_double()
#      get_vector_double()          ## (2/16/10
#      get_grid_double()
#      get_values_in_grid_double()  ## (2/17/10)
#      -----------------------------
#       Next 3 currently identical
#      -----------------------------
#      set_scalar_double()
#      set_vector_double()          ## (2/16/10)
#      set_grid_double()
#      set_values_in_grid_double()  ## (2/17/10)
#      ---------------------
#      get_scalar_long()
#      get_vector_long()
#      get_grid_long()
#      get_values_in_grid_long()    ## (2/17/10)
#      ---------------------
#      set_scalar_long()            ## (2/16/10)
#      set_vector_long()            ## (2/16/10)
#      set_grid_long()              ## (2/16/10)
#      set_values_in_grid_long()    ## (2/17/10)
#      ---------------------
#      get_input_items()
#      get_output_items()
#      ---------------------
#      get_gui_info()           ## (commented out)
#      read_gui_info_file()     ## (11/13/09)
#      get_user_input()         ## (9/24/09, 11/13/09)
#      read_config_gui()        ## (10/1/09)
#      read_config_file()       ## (5/17/10)
#      ---------------------
#      go()
#      run_model()
#      ### read_cfg_file()   # (template)
#      ### initialize        # (template)
#      ### update()          # (template)
#      finalize()
#      -------------------------
#      set_directory()                    ## obsolete soon ?
#      initialize_config_vars()           ## (5/6/10)
#      open_new_sync_file()               ## (5/7/10)
#      write_sync_file()                  ## (5/7/10)
#      read_sync_file()                   ## (5/7/10)
#      set_computed_input_vars            ## (5/6/10) over-ridden by comps.

#      -----------------------------------------------
#      These methods are not part of "IRF" interface
#      but are used by the TopoFlow components.
#      -----------------------------------------------
#      initialize_required_components()
#      get_cca_port()
#      get_cca_ports()
#      release_cca_ports()
#      ------------------------
#      add_child_port()
#      get_port_data()     # (rename to get_port_double ??)
#      set_port_data()     # (rename to set_port_double ??)
#      ------------------------
#      get_rank()
#      get_size()
#      ------------------------
#      read_grid_info()
#      initialize_basin_vars()
#      initialize_time_vars()
#      update_time()
#      print_time_and_value()
#      print_run_time()
#      print_final_report()   # (6/30/10)
#      print_traceback()  # (10/10/10)
#
#-----------------------------------------------------------------------

from numpy import *
import numpy

import sys, os, time
import traceback        # (10/10/10)

#--------------------------------------------
# (5/14/10. Can't be here because basins.py
# has import CSDMS_base at top (circular).
# See initialize_basin_vars() below.
#--------------------------------------------
# import basins
  
import cfg_files as cfg
import pixels
import rti_files
import tf_utils

#-----------------------------------------------------------------------
def unit_test():

    #-----------------------------------------
    # This function adjusts for the platform
    # and can be changed in "tf_utils.py".
    #-----------------------------------------
    # in_directory = tf_utils.TF_Test_Directory()

    c = CSDMS_component()
    c.CCA   = False
    c.DEBUG = True
    print 'Instantiated component.'

    #--------------------------------------------
    # Check the "read_gui_info_file()" function
    #--------------------------------------------
##    gui_info_file = '/Applications/Erode/gui_info/Erode_GUI_Info.txt'
##    gui_info_file = '/data/progs/erode/3.0/gui_info/Erode_GUI_Info.txt'
    gui_info_dir  = '/Users/peckhams/Desktop/topoflow3/TF_GUI_Info_Files_5_11_10/'
    gui_info_file = gui_info_dir + 'GC2D_GUI_Info.txt'

    #----------------------------------
    # Read GUI info file to get lists
    #----------------------------------
    var_names, labels, values, min_vals, max_vals, \
               desc, group_names, group_sizes = \
                   c.read_gui_info_file( gui_info_file )
    
    print 'var_names ='
    print var_names
    print 'labels ='
    print labels
    print 'values ='
    print values
    print 'min_vals ='
    print min_vals
    print 'max_vals ='
    print max_vals
    print 'desc ='
    print desc
    print 'group_names ='
    print group_names
    print 'group_sizes ='
    print group_sizes    
    return

    #-------------------------------------
    # Test the "print_run_time()" method
    #-------------------------------------
##    print ' '
##    print 'Testing "print_run_time()"...'
##    c.print_run_time(seconds=1.618)
##    c.print_run_time(seconds=60)
##    c.print_run_time(seconds=3600)
##    c.print_run_time(seconds=3600 * 24)
##    c.print_run_time(seconds=3600 * 24 * 365)
    
    #---------------------------
    # Test some of the methods
    #---------------------------
    c.a = numpy.float64(5)
    print 'c.a = numpy.float64(5)'
    print "c.is_scalar('a') =", c.is_scalar('a')
    print "c.is_grid('a')   =", c.is_grid('a')
    v1 = c.get_port_data('a', c)
    print "c.get_port_data('a',c) =", v1
    print ' '
    #-------------------------------------------------
    c.b = numpy.zeros((3,3), dtype='Float64')
    print "c.b = numpy.zeros((3,3), dtype='Float64')"
    print "c.is_scalar('b') =", c.is_scalar('b')
    print "c.is_grid('b')   =", c.is_grid('b')
    v2 = c.get_port_data('b', c)
    print "c.get_port_data('b',c) =", v2
    print ' '
    #-------------------------------------------------
    print "c.is_scalar('b[1]') =", c.is_scalar('b[1]')
    print "c.is_grid('b[1]')   =", c.is_grid('b[1]')
    v3 = c.get_port_data('b[1]', c)
    print "c.get_port_data('b[1]',c) =", v3
    print ' '
    #-------------------------------------------------    
    
    # This component has no initialize() method
##    c.initialize( mode='driver' )
##    print 'nx =', c.nx
##    print 'ny =', c.ny

#   unit_test()
#-----------------------------------------------------------------------
##def sync_file_test():
##
##    comp = CSDMS_component()
##
##    #---------------------------
##    # Sample info for a Driver
##    #---------------------------
##    comp.mode           = 'driver'
##    comp.in_directory   = os.getcwd()
##    comp.out_directory  = os.getcwd()
##    comp.site_prefix    = 'Treynor'
##    comp.case_prefix    = 'Case5'
##    
##    comp.write_sync_file( REPORT=True )
##    comp.read_sync_file(  REPORT=True )
##
###   sync_file_test()   
#-----------------------------------------------------------------------
class CSDMS_component:

    def __init__(self):

        self.CCA    = tf_utils.TF_Use_CCA()
        ## self.DEBUG     = True
        self.DEBUG       = False
        self.SKIP_ERRORS = False
        self.SILENT      = False
        self.REPORT      = False
        self.DONE        = False
        self.status      = 'created'   # (OpenMI 2.0 conventions)
        
        self.USE_GUI_SETTINGS = False
        self.in_directory     = None
        self.out_directory    = None
        self.site_prefix      = None
        self.case_prefix      = None
        self.comp_status      = 'Enabled'
        
        # NB! This probably won't work here, because a driver
        #     may be instantiated later that then changes the
        #     current working directory.
##        self.cfg_directory  = os.getcwd()
##        self.cfg_prefix     = 'Case5'
        
    #   __init__()
    #-------------------------------------------------------------------
    def get_status(self):

        #-----------------------------------------------------
        # Notes: Return component status as a string.  The
        #        possible return values are from OpenMI 2.0:
        #
        #           created, initializing, initialized,
        #           updating, updated, finalizing, finalized,
        #           failed (could add "stopped").
        #-----------------------------------------------------
        return self.status

    #   get_status()
    #-------------------------------------------------------------------
    def is_scalar(self, var_name):

        #------------------------------------------------
        # NB!  Case in var_name must be an exact match.
        #-------------------------------------------------      
        exec("n = numpy.rank(self." + var_name + ")")       
        return (n == 0)
    
    #   is_scalar()
    #-------------------------------------------------------------------
    def is_vector(self, var_name):

        #------------------------------------------------
        # NB!  Case in var_name must be an exact match.
        #------------------------------------------------     
        exec("n = numpy.rank(self." + var_name + ")")       
        return (n == 1)
    
    #   is_vector()
    #-------------------------------------------------------------------
    def is_grid(self, var_name):

        #------------------------------------------------
        # NB!  Case in var_name must be an exact match.
        #------------------------------------------------ 

        #-------------------------------------------------
        # (9/29/09) This might be causing a problem with
        # the c++ bindings for this CCA component.
        #-------------------------------------------------         
##        exec("type_str = str(type(self." + var_name + "))")
##        p1 = type_str.find("ndarray")
##        p2 = type_str.find("float")
##        if (p1 == -1) and (p2 == -1):
##            print 'ERROR: type(' + var_name + ') =' + type_str
##            return False
        #-------------------------------------------------
        # (9/29/09) This might be causing a problem with
        # the c++ bindings for this CCA component.
        #-------------------------------------------------        
##        if ("ndarray" not in type_str) and \
##           ("float" not in type_str):
##            print 'ERROR: type(' + var_name + ') =' + type_str
##            return False
        #-------------------------------------------------------        
        exec("n = numpy.rank(self." + var_name + ")")
        return (n == 2)

    #   is_grid()
    #-------------------------------------------------------------------
    def has_variable(self, var_name):

        #------------------------------------------------------
        # If var_name includes square brackets for subscripts
        # remove them to get the actual variable name.
        #------------------------------------------------------
        bracket_pos = var_name.find('[')
        if (bracket_pos != -1):
            key = var_name[0:bracket_pos]
        else:
            key = var_name

        #---------------------------------------------------
        # Does current component have requested variable ?
        #---------------------------------------------------
        VARIABLE_FOUND = self.__dict__.has_key(key)
        if not(VARIABLE_FOUND):
            print 'ERROR: Component does not have the'
            print '       requested variable: ' + var_name
            print ' '
            
        return VARIABLE_FOUND
        
    #   has_variable()
    #-------------------------------------------------------------------        
    def get_scalar_double(self, var_name):

        #------------------------------------
        # Note: The next line doesn't work.
        #------------------------------------
        ## exec("return self." + var_name)

        #---------------------------------------------------
        # Does current component have requested variable ?
        #---------------------------------------------------
        # This is not used yet because possible impact on
        # performance has not be tested yet. (2/17/10)
        # If it does get used later, it will be added to
        # all of the "getters".
        #---------------------------------------------------        
##        if not(self.has_variable(var_name)):
##            return float64(0)
 
        try:
            exec("result = self." + var_name)
            return numpy.float64(result)
        except:
            print 'ERROR in CSDMS_base.get_scalar_double().'
            print '    Returning 0.'
            return numpy.float64(0)

            ############## flush output here ??
        
    #   get_scalar_double()
    #-------------------------------------------------------------------
    def get_vector_double(self, var_name):

        #---------------------------------------------------------
        # Note: This function was causing a "segmentation fault
        #       in gui-backend.sh" error message when trying to
        #       run TopoFlow through the CMT (in CCA framework).
        #       Solution was to use numpy.array, as shown.
        #       (2/17/10)
        #---------------------------------------------------------
        try:
            exec("result = self." + var_name)
            return numpy.array(result, dtype='float64')
            #-------------------------
            # NB! This doesn't work.
            #-------------------------
            # return numpy.float64(result)
        except:
            print 'ERROR in CSDMS_base.get_vector_double().'
            print '    Returning zeros.'
            return numpy.zeros([1], dtype='float64')
        
    #   get_vector_double()
    #-------------------------------------------------------------------
    def get_grid_double(self, var_name):

        try:
            exec("result = self." + var_name)
            return numpy.float64(result)
        except:
            print 'ERROR in CSDMS_base.get_grid_double().'
            print '    Returning zeros.'
            return numpy.zeros([1,1], dtype='float64')
        
    #   get_grid_double()
    #-------------------------------------------------------------------
    def get_values_in_grid_double(self, var_name, IDs):

        #---------------------------------------------------------
        # Note: This function was causing a "segmentation fault
        #       in gui-backend.sh" error message when trying to
        #       run TopoFlow through the CMT (in CCA framework).
        #       Solution was to use numpy.array, as shown.
        #       (2/18/10)
        #---------------------------------------------------------
        # Notes: This function was tested in the new Diversions
        #        component on (2/18/10).
        #---------------------------------------------------------
        try:
            exec("result = self." + var_name + '.flat[IDs]')
            return numpy.array(result, dtype='float64')
            ## return numpy.float64(result)
        except:
            print 'ERROR in CSDMS_base.get_values_in_grid_double().'
            print '    Returning zeros.'
            return numpy.zeros(len(IDs), dtype='float64')
        
    #   get_values_in_grid_double()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def set_scalar_double(self, var_name, scalar):

        exec("self." + var_name + " = numpy.float64(scalar)")
    
    #   set_scalar_double()
    #-------------------------------------------------------------------
    def set_vector_double(self, var_name, vector):

        #--------------------------------------------------
        # Notes: First method here should be more robust.
        #        See Notes for get_vector_double().
        #--------------------------------------------------
        exec("self." + var_name + " = numpy.array(vector, dtype='float64')")

        #-----------------------------------
        # Original method (before 2/17/10)
        #-----------------------------------        
        # exec("self." + var_name + " = numpy.float64(vector)")
    
    #   set_vector_double()
    #-------------------------------------------------------------------
    def set_grid_double(self, var_name, grid):
         
        exec("self." + var_name + " = numpy.float64(grid)")
    
    #   set_grid_double()
    #-------------------------------------------------------------------
    def set_values_in_grid_double(self, var_name, IDs, values):

        #--------------------------------------------------------
        # Notes: This function was tested in the new Diversions
        #        component on (2/18/10).
        #--------------------------------------------------------

        exec("self." + var_name + ".flat[IDs] = values")
        
        # exec("self." + var_name + ".flat[IDs] = numpy.float64(values)")
        
    #   set_values_in_grid_double()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def get_scalar_long(self, var_name):

        exec("result = numpy.int32(self." + var_name + ")")
        return result
    
    #   get_scalar_long()
    #-------------------------------------------------------------------
    def get_vector_long(self, var_name):

        #--------------------------------------------
        # Notes: See Notes for get_vector_double().
        #--------------------------------------------
        try:
            exec("result = self." + var_name)
            return numpy.array(result, dtype='int32')
            #-------------------------
            # NB! This doesn't work.
            #-------------------------
            # return numpy.int32(result)
        except:
            print 'ERROR in CSDMS_base.get_vector_long().'
            print '    Returning zeros.'
            return numpy.zeros([1], dtype='int32')
        
    #   get_vector_long()
    #-------------------------------------------------------------------
    def get_grid_long(self, var_name):

        exec("result = numpy.int32(self." + var_name + ")")
        return result
    
    #   get_grid_long()
    #-------------------------------------------------------------------
    def get_values_in_grid_long(self, var_name, IDs):

        try:
            exec("result = self." + var_name + '.flat[IDs]')
            return numpy.int32(result)
        except:
            print 'ERROR in CSDMS_base.get_values_in_grid_long().'
            print '    Returning zeros.'
            return numpy.zeros(len(IDs), dtype='int32')
        
    #   get_values_in_grid_long()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def set_scalar_long(self, var_name, scalar):

        exec("self." + var_name + " = numpy.int32(scalar)")
    
    #   set_scalar_long()
    #-------------------------------------------------------------------
    def set_vector_long(self, var_name, vector):

        #--------------------------------------------------
        # Notes: First method here should be more robust.
        #        See Notes for get_vector_double().
        #--------------------------------------------------
        exec("self." + var_name + " = numpy.array(vector, dtype='int32')")

        #-----------------------------------
        # Original method (before 2/17/10)
        #-----------------------------------
        # exec("self." + var_name + " = numpy.int32(vector)")
        
    #   set_vector_long()
    #-------------------------------------------------------------------
    def set_grid_long(self, var_name, grid):

        exec("self." + var_name + " = numpy.int32(grid)")
    
    #   set_grid_long()
    #-------------------------------------------------------------------
    def set_values_in_grid_long(self, var_name, IDs, values):

        #----------------------------------------------------------
        # Note: Type of "values" should already be long (SIDL).
        #----------------------------------------------------------
        # exec("self." + var_name + ".flat[IDs] = numpy.int32(values)")
        exec("self." + var_name + ".flat[IDs] = values")
        
    #   set_values_in_grid_long()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def get_input_items(self):

        #-------------------------------------------------------------
        # Note: There may be a way to retrieve this list auto-
        #       matically using code similar to self.has_variable().
        #-------------------------------------------------------------
        items = ['None']
        return numpy.array(items)   # (string array vs. list)
    
    #   get_input_items()
    #-------------------------------------------------------------------
    def get_output_items(self):

        items = ['None']
        return numpy.array(items)   # (string array vs. list)
    
    #   get_output_items()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    def read_gui_info_file(self, gui_info_file):

        #------------------------------------------
        # Open file to read and skip header lines
        #------------------------------------------
        file_unit = open(gui_info_file, 'r')
        ## line1 = file_unit.readline()
        cfg.skip_header( file_unit, n_lines=4 )

        group_names = []
        group_sizes = []
        group_size  = 0
        #------------------
        var_names   = []
        labels      = []
        types       = []
        values      = []
        min_vals    = []
        max_vals    = []
        desc        = []

        #-----------------------------
        # Read information from file
        #-----------------------------
        while (True):
            line = file_unit.readline()
            if (line == ''):
                break
            words = line.split('|')
            char1 = words[0][0]  # (don't use .strip() due to "\n")
            if (len(words) > 6) and (char1 != '#'):
                group_size += 1  #####
                #--------------------------------------
                var_names.append( words[0].strip() )
                labels.append( words[1].strip() )
                #--------------------------------------
                vtype   = words[2].strip().lower()
                value   = words[3].strip()
                min_val = words[4].strip()
                max_val = words[5].strip()
                if (vtype != 'string'):
                    value   = eval( value )
                    min_val = eval( min_val )
                    max_val = eval( max_val )
                if (vtype in ['int', 'long']):
                    #-------------------------------
                    # Need this to avoid a warning
                    #-------------------------------
                    value   = numpy.int32( value )
                    min_val = numpy.int32( min_val )
                    max_val = numpy.int32( max_val )
                types.append( vtype )
                values.append( value )
                min_vals.append( min_val )
                max_vals.append( max_val )
                #--------------------------------------
                desc.append( words[6].strip() )
            elif (len(words) == 2) and \
                 (words[0].strip().upper() == 'PPF_GROUP_NAME'):
                group_names.append( words[1].strip() )
                if (group_size > 0):
                    group_sizes.append( group_size )
                    group_size = 0
        group_sizes.append( group_size ) # (last one)

        #--------------
        # For testing
        #--------------
##        print 'group_names =', group_names
##        print 'group_sizes =', group_sizes
##        print ' '
##        print 'var_names =', var_names
##        print 'labels    =', labels
##        print 'types     =', types
##        print 'values    =', values
##        print 'min_vals  =', min_vals
##        print 'max_vals  =', max_vals
            
        file_unit.close()
    
        return var_names, labels, types, values, min_vals, max_vals, \
               desc, group_names, group_sizes

    #   read_gui_info_file()
    #-------------------------------------------------------------------
    def get_gui_info(self):

        #-----------------------------------------------------
        # (6/22/10) Used for basins.py in read_config_file().
        #-----------------------------------------------------
        # Notes: Each component should over-ride this method
        #        with one that provides its gui info.
        #-----------------------------------------------------
        directory = "/data/progs/topoflow/3.1/gui_info/"
        file_name = "Null_Component.cfg"
        self.gui_info_file = (directory + file_name)
        # self.cfg_file = (directory + file_name)
        self.dialog_title  = "Component Parameters"

    #   get_gui_info()
    #-------------------------------------------------------------------
    def get_user_input( self, services, d_services,
                        button_label="Configure",
                        dialog_title="Component Parameters",
                        gui_info_file=None ):

        #-------------------------------------------------------
        # Note: This uses a CCA "Parameter Port" to launch a
        #       GUI dialog that prompts for input parameters
        #       to the current component instance.

        #       A call to this method can be inserted into
        #       the setServices() method of a CCA component's
        #       impl file.

        #       New button_label default: "Configure" (5/6/10)
        #-------------------------------------------------------
        import gov.cca.ports.ParameterPortFactory

        #--------------------------------------
        # Get info about this component's GUI
        #--------------------------------------
        try:
            #-----------------------------------------------
            # Ignore gui_info_file and dialog_title args
            # which will be removed from Impl files later.
            #-----------------------------------------------
            self.get_gui_info()
            gui_info_file = self.gui_info_file
            dialog_title  = self.dialog_title
            print 'Reading GUI info from: ' + gui_info_file
        except:
            print 'ERROR in CSDMS_base.get_user_input():'
            print '    gui_info_file not found.'
            return
        
        #----------------------------------
        # Read lists from "gui_info_file"
        #----------------------------------        
        var_names, labels, types, \
            values, min_vals, max_vals, \
            desc, group_names, group_sizes = \
            self.read_gui_info_file( gui_info_file )
        
        #----------------------------------------------
        # self.userinput = d_services.createTypeMap() 
        # Do we really need to store this in self ?
        # Seems it is only a local variable.
        # (self.userinput -> dialog)
        #----------------------------------------------
        dialog = d_services.createTypeMap()    

        try:
            port = d_services.getPort("ppf")
        except sidl.BaseException._Exception, e:
            port = None

        if (self.DEBUG):
            if (port == None):
                print 'FAILURE: Unable to get generic CCA port'
            else:
                print 'SUCCESS: Got a generic CCA port.'
        
        if not port:
            print 'getPort("ppf") returned nil'
        else: 
            ppf = gov.cca.ports.ParameterPortFactory.ParameterPortFactory( port ) 
            if not ppf:
                print 'Error casting to gov.cca.ports.ParameterPortFactory'
            else:
                if (self.DEBUG): print 'SUCCESS: Cast port to PPF port.'
            
        ppf.initParameterData(dialog, button_label)
        ppf.setBatchTitle(dialog, dialog_title)

        if (self.DEBUG):
            print 'Starting for loop over parameter list...'
            
        #----------------------------------------------
        # Add rows to prompt for requested parameters
        #--------------------------------------------------------
        # Note that var_names must be unique or they won't show
        #--------------------------------------------------------
        group_num = 0
        group_sum = 0
        for k in xrange(len(var_names)):
##            if (self.DEBUG):
##                print 'k =', k   #############
##                print 'types[k].lower() =', types[k].lower()
            
            #----------------------------------------------
            # Start a new "panel" with new "group name" ?
            #----------------------------------------------
            if (k == group_sum):
                ppf.setGroupName( dialog, group_names[ group_num ] )
                group_sum += group_sizes[ group_num ]
                group_num += 1
##                print 'group_num =', group_num
##                print 'group_sum =', group_sum
                
            #------------------------------------------
            # Prompt for a variable of requested type
            #------------------------------------------
            type_k = types[k].lower()
            if   (type_k in ['double', 'float']):
                ppf.addRequestDouble(dialog, var_names[k],
                                     desc[k], labels[k], values[k],
                                     min_vals[k], max_vals[k])
            elif (type_k in ['int', 'long']):
                ppf.addRequestInt(dialog, var_names[k],
                                  desc[k], labels[k], values[k],
                                  min_vals[k], max_vals[k])
            elif (type_k == 'string'):
                ppf.addRequestString(dialog, var_names[k],
                                     desc[k], labels[k], values[k] )
                                     ### min_vals[k], max_vals[k])
            else:
                print 'ERROR in CSDMS_base.get_user_input().'
                print '   Unsupported type name: ' + type_k
                
        if (self.DEBUG):
            print 'Finished with for loop.'
            print 'Calling ppf.addParameterPort()...'
        
        ppf.addParameterPort(dialog, services)
        d_services.releasePort("ppf")  
        if (self.DEBUG):
            #-------------------------------------------------
            # Note: This message is printed even if the user
            #       has not yet opened the PPF dialog.
            #-------------------------------------------------
            print 'Finished with get_user_input().'
        
        #-----------------------------------------------
        # Does call to "releasePort()" imply a closing
        # of the dialog ???
        #-----------------------------------------------
        # Maybe need to sleep for a bit here.
        # time.sleep(0.2)
        
        #---------------------------------------------------------
        # Try to store the user-entered values directly into
        # the current component's state.  This probably requires
        # that user has dismissed dialog with "OK" button.
        #---------------------------------------------------------
        # If we can't do this here, then wrap this block into a
        # small function that can be called by "go()" method.
        #---------------------------------------------------------
        self.dialog = dialog
        ## self.read_config_gui(var_names, types, dialog_str='dialog')

        #--------------------------------------------------------
        # Store these for use by "CSDMS_base.read_config_gui()"
        #--------------------------------------------------------
        # Note: This function is called (and DEBUG message is
        #       printed, see above) even if user has not
        #       opened the PPF dialog.  This is the only place
        #       where USE_GUI_SETTINGS is set to True, and this
        #       function is only called from CCA Impl files.
        #----------------------------------------------------------  
        self.PPF_var_names = var_names
        self.PPF_types     = types
        self.USE_GUI_SETTINGS = True
        
    #   get_user_input()
    #-------------------------------------------------------------------
    def read_config_gui(self, var_names=None, types=None,
                        dialog_str='self.dialog'):

        if not(self.SILENT):
            print 'Saving GUI settings into component state.'

        if (var_names == None):
            var_names = self.PPF_var_names
        if (types == None):
            types = self.PPF_types

        last_var_name = ''
        
        #-----------------------------------------
        # Save user input into component's state
        #-----------------------------------------
        for k in xrange(len(var_names)):
            var_name      = var_names[k]
            var_type      = types[k].lower()
            READ_SCALAR   = False
            READ_FILENAME = False

            #----------------------------------------------
            # Does var_name end with an array subscript ?
            #----------------------------------------------
            p1 = var_name.rfind('[')
            p2 = var_name.rfind(']')
            if (p1 > 0) and (p2 > p1):
                var_base  = var_name[:p1]
                subscript = var_name[p1:p2+1]
                var_name_file_str = var_base + '_file' + subscript
            else:
                var_base = var_name
                var_name_file_str = var_name + '_file'

            #--------------------------------------------
            # Update var_type based on droplist setting
            #--------------------------------------------
            if (last_var_name.startswith(var_base + '_type')):
                exec( "type_choice = self." + last_var_name )
                if (type_choice.lower() == 'scalar'):
                    #--------------------------------------------------
                    # It seems that things will work as long as the
                    # "type" and "value" fields in the GUI CFG file
                    # are consistent.  Don't change var_type here.
                    #
                    # Otherwise get this message:
                    # "Mismatch with value found in typemap
                    #  (requested type String, actual type Double)."
                    #--------------------------------------------------
                    exec( "self." + var_name_file_str + " = ''")
                    READ_SCALAR = True
                    ## var_type = 'double'
                else:
                    exec( "self." + var_name + " = 0.0")
                    READ_FILENAME = True
                    ## var_type = 'string'
                    
            # print 'var_name, var_type =', var_name, ',', var_type
        
            #-----------------------------------           
            # Read a value of type "var_type"
            #-----------------------------------
            # Convert scalars to numpy scalars
            #-----------------------------------
            if (var_type in ['double', 'float']):
                exec( "value = " + dialog_str + ".getDouble( var_name, 0.0 )" )
                value = numpy.float64( value )
                exec( "self." + var_name + " = value" )
            elif (var_type in ['long', 'int']):
                exec( "value = " + dialog_str + ".getInt( var_name, 0 )" )
                value = numpy.int32( value )
                exec( "self." + var_name + " = value" )
            elif (var_type == 'string'):
                #-----------------------------------------
                # Get the value string for this var_name
                #----------------------------------------------------
                # If string has a placeholder filename prefix, then
                # expand it here.  Need to use original "var_name"
                # without appending "_file" until assignment.
                #----------------------------------------------------
                exec("s = " + dialog_str + ".getString( var_name, 'ERROR' )" )
                if (s[:13] == '<case_prefix>'):
                    value_str = (self.case_prefix + s[13:])
                elif (s[:13] == '<site_prefix>'):
                    value_str = (self.site_prefix  + s[13:])
                else:
                    value_str = s

                #-----------------------------------------------
                # If var_name starts with "SAVE_" and value is
                # Yes or No, then convert to Python boolean.
                #-----------------------------------------------
                if (var_name[:5] == 'SAVE_'):
                    VALUE_SET = True
                    if (s.lower() == 'yes'):
                        exec( "self." + var_name + " = True" )
                    elif (s.lower() == 'no'):
                        exec( "self." + var_name + " = False" )
                    else:
                        VALUE_SET = False
                else:
                    VALUE_SET = False
                #----------------------------------------------------------
                if not(VALUE_SET):
                    if (READ_FILENAME):
                        exec( "self." + var_name_file_str + " = value_str" )
                    elif (READ_SCALAR):
                        exec( "self." + var_name + " = numpy.float64(value_str)")
                    else:
                        exec( "self." + var_name + " = value_str" )
                    
            else:
                print 'ERROR in CSDMS_base.read_config_gui().'
                print '   Unsupported data type = ' + var_type + '.'
                print ' '

            last_var_name = var_name
            
    #   read_config_gui()
    #-------------------------------------------------------------------
    def read_config_file(self):

        if not(self.SILENT):
            print 'Reading config file into component state.'

        #---------------------------
        # Get name of the cfg_file
        #---------------------------
        self.get_gui_info()
        cfg_extension = self.get_cfg_extension()
        cfg_directory = (os.getcwd() + os.sep)
        file_name     = (self.cfg_prefix + cfg_extension)
        self.cfg_file = (cfg_directory + file_name)
        if (self.DEBUG):
            print 'cfg_file =', self.cfg_file
        if not(os.path.exists(self.cfg_file)):
            print 'WARNING: cfg_file not found:'
            print '         ' + self.cfg_file
            return
         
        #----------------------------------
        # Read lists from "gui_info_file"
        #----------------------------------
        gui_info_file = self.cfg_file   #######
        if (gui_info_file == None):
            print 'ERROR in CSDMS_base.read_config_file():'
            print '    gui_info_file not found.'
            return
        var_names, labels, types, \
            values, min_vals, max_vals, \
            desc, group_names, group_sizes = \
            self.read_gui_info_file( gui_info_file )

        last_var_name = ''
        
        #-----------------------------------------
        # Save user input into component's state
        #-----------------------------------------
        for k in xrange(len(var_names)):
            var_name      = var_names[k]
            var_type      = types[k].lower()
            value         = values[k]  # (list of diff. types)
            READ_SCALAR   = False
            READ_FILENAME = False

            #----------------------------------------------
            # Does var_name end with an array subscript ?
            #----------------------------------------------
            p1 = var_name.rfind('[')
            p2 = var_name.rfind(']')
            if (p1 > 0) and (p2 > p1):
                var_base  = var_name[:p1]
                subscript = var_name[p1:p2+1]
                var_name_file_str = var_base + '_file' + subscript
            else:
                var_base = var_name
                var_name_file_str = var_name + '_file'

            #--------------------------------------------
            # Update var_type based on droplist setting
            #--------------------------------------------
            if (last_var_name.startswith(var_base + '_type')):
                exec( "type_choice = self." + last_var_name )
                if (type_choice.lower() == 'scalar'):
                    #--------------------------------------------------
                    # It seems that things will work as long as the
                    # "type" and "value" fields in the GUI CFG file
                    # are consistent.  Don't change var_type here.
                    #
                    # Otherwise get this message:
                    # "Mismatch with value found in typemap
                    #  (requested type String, actual type Double)."
                    #--------------------------------------------------
                    exec( "self." + var_name_file_str + " = ''")
                    READ_SCALAR = True
                    ## var_type = 'double'
                else:
                    exec( "self." + var_name + " = 0.0")
                    READ_FILENAME = True
                    ## var_type = 'string'

            # print 'var_name, var_type =', var_name, ',', var_type

            #-----------------------------------           
            # Read a value of type "var_type"
            #-----------------------------------
            # Convert scalars to numpy scalars
            #-----------------------------------
            if (var_type in ['double', 'float']):
                value = numpy.float64( value )
                exec( "self." + var_name + " = value" )
            elif (var_type in ['long', 'int']):
                value = numpy.int32( value )
                exec( "self." + var_name + " = value" )
            elif (var_type == 'string'):
                #-----------------------------------------
                # Get the value string for this var_name
                #----------------------------------------------------
                # If string has a placeholder filename prefix, then
                # expand it here.  Need to use original "var_name"
                # without appending "_file" until assignment.
                #----------------------------------------------------
                s = value
                if (s[:13] == '<case_prefix>'):
                    value_str = (self.case_prefix + s[13:])
                elif (s[:13] == '<site_prefix>'):
                    value_str = (self.site_prefix  + s[13:])
                else:
                    value_str = s

                #-----------------------------------------------
                # If var_name starts with "SAVE_" and value is
                # Yes or No, then convert to Python boolean.
                #-----------------------------------------------
                if (var_name[:5] == 'SAVE_'):
                    VALUE_SET = True
                    if (s.lower() == 'yes'):
                        exec( "self." + var_name + " = True" )
                    elif (s.lower() == 'no'):
                        exec( "self." + var_name + " = False" )
                    else:
                        VALUE_SET = False
                else:
                    VALUE_SET = False
                #----------------------------------------------------------
                if not(VALUE_SET):
                    if (READ_FILENAME):
                        exec( "self." + var_name_file_str + " = value_str" )
                    elif (READ_SCALAR):
                        exec( "self." + var_name + " = numpy.float64(value_str)")
                    else:
                        exec( "self." + var_name + " = value_str" )
            else:
                print 'ERROR in CSDMS_base.read_config_file().'
                print '   Unsupported data type = ' + var_type + '.'
                print ' '

            last_var_name = var_name

    #   read_config_file()
    #-------------------------------------------------------------------
    def read_config_file2(self):

        #----------------------------------------------------------
        # Notes: This version reads configuration settings from
        #        a new type of CFG file that only has var_name,
        #        value and type; more like key-value. (5/9/11)
        #        The GUI is no longer generated using a CFG file.
        #----------------------------------------------------------
        
        if not(self.SILENT):
            print 'Reading config file into component state.'

        #---------------------------
        # Get name of the cfg_file
        #---------------------------
        ### self.get_gui_info()
        cfg_extension = self.get_cfg_extension()
        cfg_directory = (os.getcwd() + os.sep)
        file_name     = (self.cfg_prefix + cfg_extension)
        self.cfg_file = (cfg_directory + file_name)
        if (self.DEBUG):
            print 'cfg_file =', self.cfg_file
        if not(os.path.exists(self.cfg_file)):
            print 'WARNING: cfg_file not found:'
            print '         ' + self.cfg_file
            return

        #-----------------------------
        # Open CFG file to read data
        #-----------------------------
        cfg_unit = open( self.cfg_file, 'r' )
        last_var_name = ''

        #-----------------------------------------
        # Save user input into component's state
        #--------------------------------------------------
        # Recall that a "blank line", with just a (hidden)
        # newline character will not be null and will
        # have len(line) = 1.
        #--------------------------------------------------
        while (True):
            line  = cfg_unit.readline()
            if (line == ''):
                break                  # (reached end of file)
            COMMENT = (line[0] == '#')
            #--------------------------------------------
            # Using "|" as a delimiter means we can use
            # " ", "," or "=" in the filenames.
            #--------------------------------------------
            words   = line.split('|')  # (split on equals)
            if (len(words) == 3) and not(COMMENT):
                var_name = words[0]
                var_type = words[1]
                value    = words[2]
                READ_SCALAR   = False
                READ_FILENAME = False

                #----------------------------------------------
                # Does var_name end with an array subscript ?
                #----------------------------------------------
                p1 = var_name.rfind('[')
                p2 = var_name.rfind(']')
                if (p1 > 0) and (p2 > p1):
                    var_base  = var_name[:p1]
                    subscript = var_name[p1:p2+1]
                    var_name_file_str = var_base + '_file' + subscript
                else:
                    var_base = var_name
                    var_name_file_str = var_name + '_file'

                #--------------------------------------------
                # Update var_type based on droplist setting
                #--------------------------------------------
                if (last_var_name.startswith(var_base + '_type')):
                    exec( "type_choice = self." + last_var_name )
                    if (type_choice.lower() == 'scalar'):
                        #--------------------------------------------------
                        # It seems that things will work as long as the
                        # "type" and "value" fields in the GUI CFG file
                        # are consistent.  Don't change var_type here.
                        #
                        # Otherwise get this message:
                        # "Mismatch with value found in typemap
                        #  (requested type String, actual type Double)."
                        #--------------------------------------------------
                        exec( "self." + var_name_file_str + " = ''")
                        READ_SCALAR = True
                        ## var_type = 'double'
                    else:
                        exec( "self." + var_name + " = 0.0")
                        READ_FILENAME = True
                        ## var_type = 'string'

                # print 'var_name, var_type =', var_name, ',', var_type

                #-----------------------------------           
                # Read a value of type "var_type"
                #-----------------------------------
                # Convert scalars to numpy scalars
                #-----------------------------------
                if (var_type in ['double', 'float']):
                    value = numpy.float64( value )
                    exec( "self." + var_name + " = value" )
                elif (var_type in ['long', 'int']):
                    value = numpy.int32( value )
                    exec( "self." + var_name + " = value" )
                elif (var_type == 'string'):
                    #-----------------------------------------
                    # Get the value string for this var_name
                    #----------------------------------------------------
                    # If string has a placeholder filename prefix, then
                    # expand it here.  Need to use original "var_name"
                    # without appending "_file" until assignment.
                    #----------------------------------------------------
                    s = value
                    if (s[:13] == '<case_prefix>'):
                        value_str = (self.case_prefix + s[13:])
                    elif (s[:13] == '<site_prefix>'):
                        value_str = (self.site_prefix  + s[13:])
                    else:
                        value_str = s

                    #-----------------------------------------------
                    # If var_name starts with "SAVE_" and value is
                    # Yes or No, then convert to Python boolean.
                    #-----------------------------------------------
                    if (var_name[:5] == 'SAVE_'):
                        VALUE_SET = True
                        if (s.lower() == 'yes'):
                            exec( "self." + var_name + " = True" )
                        elif (s.lower() == 'no'):
                            exec( "self." + var_name + " = False" )
                        else:
                            VALUE_SET = False
                    else:
                        VALUE_SET = False
                    #----------------------------------------------------------
                    if not(VALUE_SET):
                        if (READ_FILENAME):
                            exec( "self." + var_name_file_str + " = value_str" )
                        elif (READ_SCALAR):
                            exec( "self." + var_name + " = numpy.float64(value_str)")
                        else:
                            exec( "self." + var_name + " = value_str" )
                else:
                    print 'ERROR in CSDMS_base.read_config_file().'
                    print '   Unsupported data type = ' + var_type + '.'
                    print ' '

                last_var_name = var_name

    #   read_config_file2() 
    #-------------------------------------------------------------------
##    def read_config_gui_last(self, var_names=None, types=None,
##                             dialog_str='self.dialog'):
##
##        if (self.DEBUG):
##            print 'Saving user input into component state.'
##
##        if (var_names == None):
##            var_names = self.PPF_var_names
##        if (types == None):
##            types = self.PPF_types
##
##        #-----------------------------------------
##        # Save user input into component's state
##        #-----------------------------------------
##        for k in xrange(len(var_names)):
##            if (types[k].lower() == 'double'):            
##                exec( "self." + var_names[k] + " = " +
##                      dialog_str + ".getDouble( var_names[k], 0.0 )" )
##            elif (types[k].lower() == 'float'):            
##                exec( "self." + var_names[k] + " = " +
##                      dialog_str + ".getDouble( var_names[k], 0.0 )" )
##            elif (types[k].lower() == 'long'):                
##                exec( "self." + var_names[k] + " = " +
##                      dialog_str + ".getInt( var_names[k], 0 )" )
##            elif (types[k].lower() == 'int'):                
##                exec( "self." + var_names[k] + " = " +
##                      dialog_str + ".getInt( var_names[k], 0 )" )
##            elif (types[k].lower() == 'string'):
##                #----------------------------------------
##                # If string has a placeholder filename
##                # prefix, then expand it here.
##                #----------------------------------------
##                s = var_names[k]
##                if (s[:13] == '<case_prefix>'):
##                    var_names[k] = (self.case_prefix + s[13:])
##                if (s[:13] == '<site_prefix>'):
##                    var_names[k] = (self.site_prefix  + s[13:])
##                #--------------------------------------------------
##                exec( "self." + var_names[k] + " = " +
##                      dialog_str + ".getString( var_names[k], '0.0' )" )
##            else:
##                print 'ERROR in CSDMS_base.read_config_gui().'
##                print '   Unsupported data type = ' + types[k] + '.'
##                print ' '
##    
##    #   read_config_gui()   
    #-------------------------------------------------------------------
    def go(self):

        self.run_model()

    #   go()
    #-------------------------------------------------------------------
    def run_model(self, cfg_directory=None, cfg_prefix=None,
                  n_steps=5):

        #--------------------------------------------------
        # All components, including this one (the driver)
        # will look in the CWD for their CFG file.
        #--------------------------------------------------
        if (cfg_directory != None):
            os.chdir( cfg_directory )
        self.cfg_prefix = cfg_prefix

        #-----------------------------------------
        # Initialize the model run (driver mode)
        #---------------------------------------------
        # This will set in_directory, out_directory,
        # site_prefix and case_prefix
        #---------------------------------------------
        self.initialize( cfg_prefix=cfg_prefix, mode='driver' )
            
        #----------------------------------------------------------- 
        # Note:  If the number of timesteps is specified in a
        #        component's CFG file and is then saved by
        #        "read_cfg_file()" as "n_steps" then we should
        #        honor that setting.  Otherwise we use the n_steps
        #        argument.
        #-----------------------------------------------------------    
        if (hasattr(self, 'n_steps')):
            ## print 'NUMBER OF STEPS =', self.n_steps  ####
            n_steps = self.n_steps

        #-------------------------------------------        
        # Note: __init__() sets self.DONE to False
        #-------------------------------------------
        while not(self.DONE):
            if not(self.SKIP_ERRORS):
                #-------------------------------------------
                # Exceptions will not be caught and
                # Python error messages will be displayed.
                #-------------------------------------------
                if (self.DEBUG):
                    print 'time_index =', self.time_index
                self.update()
            else:   
                try:
                    self.update()
                except:
                    print 'ERROR in run_model() method at:'
                    print '   time_index =', self.time_index
                    self.status = 'failed'
                    self.DONE = True

            #------------------------------------------------
            # If the model has set self.DONE = True, then
            # stop, even if we haven't reached n_steps yet.
            #------------------------------------------------
            TIMES_UP  = (self.time_index >= n_steps)
            self.DONE = (self.DONE or TIMES_UP)

        #-------------------------
        # Finalize the model run
        #-------------------------
        self.finalize()

    #   run_model()
    #-------------------------------------------------------------
    def initialize(self, cfg_prefix=None, mode="nondriver"):

        self.status     = 'initializing'  # (OpenMI 2.0 convention)
        self.mode       = mode
        self.cfg_prefix = cfg_prefix
        
        #-----------------------------------------------
        # Load component parameters from a config file
        # Will use cfg_file from above.
        #-----------------------------------------------
        ## self.set_constants()
        self.initialize_config_vars() 
        self.read_grid_info()
        self.initialize_time_vars()

        #---------------------------------------------
        # Open input files needed to initialize vars 
        #---------------------------------------------
        self.open_input_files()
        self.read_input_files()

##        self.open_output_files()

        self.status = 'initialized'  # (OpenMI 2.0 convention)
        
    #   initialize()
    #-------------------------------------------------------------------
##    def update(self, time_seconds=None):
##
##        self.status = 'updating'  # (OpenMI 2.0 convention)
##
##        #-----------------------------------
##        # Update computed values here with
##        # a series of "update_var()" calls
##        #----------------------------------
##        # self.update_var1()
##        # self.update_var2()
##
##        #------------------------
##        # Update internal clock
##        #------------------------
##        self.update_time()
##
##        #-------------------------------
##        # Check for NaNs, etc. in var1
##        #-------------------------------    
##        # self.check_var1()
##
##        #------------------------------------------
##        # Read next infil vars from input files ?
##        #------------------------------------------
##        self.read_input_files()
##
##        #----------------------------------------------
##        # Write user-specified data to output files ?
##        #----------------------------------------------
##        self.write_output_files(time_seconds)
##        self.status = 'updated'  # (OpenMI 2.0 convention)
##        
##    #   update()
    #-------------------------------------------------------------------
    def finalize(self):

        self.status = 'finalizing'  # (OpenMI 2.0 convention)
        self.close_input_files()    ##  TopoFlow input "data streams"
        self.close_output_files()
        self.status = 'finalized'  # (OpenMI 2.0 convention)
        
    #   finalize()
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
    #-------------------------------------------------------------------
##    def set_directory(self, in_directory, site_prefix, case_prefix):
##
##        #-------------------------------------
##        # Make sure these are always defined
##        #-------------------------------------
##        if (in_directory is None):
##            in_directory = tf_utils.Current_Directory()
##        if (site_prefix is None):
##            site_prefix = 'Null'
##            ## Get_Site_Prefix has filepath arg and
##            ## a different purpose.
##            ## site_prefix = tf_utils.Get_Site_Prefix()
##        if (case_prefix is None):
##            case_prefix = tf_utils.Get_Case_Prefix()
##
##        #--------------------------------------------------
##        # Add trailing separator to directory, if missing
##        #--------------------------------------------------
##        if (directory[-1] != os.sep):
##            directory += os.sep
##        
##        self.directory   = directory
##        self.site_prefix = site_prefix
##        self.case_prefix = case_prefix
##
##        #----------------------------
##        # CD to working directory ?
##        # May contain blank spaces
##        #----------------------------
##        os.chdir( self.in_directory )
##        
##    #   set_directory()
    #-------------------------------------------------------------------
    def check_directories(self):

        #----------------------------------------
        # Note:  Defaults are set in __init_().
        #----------------------------------------
        if (self.in_directory == None):
            self.in_directory = os.getcwd() + os.sep
            ## self.in_directory = self.cfg_directory
        #-----------------------------------------------
        if (self.out_directory == None):
            self.out_directory = os.getcwd() + os.sep
        #-----------------------------------------------
        if (self.site_prefix == None):
            self.site_prefix = self.cfg_prefix
        #-----------------------------------------------
        if (self.case_prefix == None):
            self.case_prefix = self.cfg_prefix

##        print 'self.in_directory  =', self.in_directory
##        print 'self.out_directory =', self.out_directory
        
        #--------------------------------------------------
        # Add trailing separator to directory, if missing
        #--------------------------------------------------
        # (self.in_directory != ''):
        if (self.in_directory[-1] != os.sep):
            self.in_directory += os.sep
        #----------------------------------------
        # (self.out_directory != ''):
        if (self.out_directory[-1] != os.sep):
            self.out_directory += os.sep

        #-----------------------------------------------
        # Expand path abbreviations like "~" (5/19/10)
        #-----------------------------------------------
        self.in_directory  = os.path.expanduser( self.in_directory  )
        self.out_directory = os.path.expanduser( self.out_directory )
        
    #   check_directories() 
    #-------------------------------------------------------------------
    def initialize_config_vars(self):
   
        #------------------------------------------------------
        # Notes:  In the CMT tool, the component acting as
        #         the driver is the one with a "Run" button.
        #
        #         When a user clicks on a driver's "Run"
        #         button, the CCA "go()" method calls the
        #         driver's "run_()" method.  This sets
        #         the "mode" to "driver" and passes the
        #         mode to initialize().  If the initialize()
        #         method is called by another component (e.g.
        #         a driver), then the mode keeps the default
        #         setting of "nondriver".
        #
        #         Only the driver can change the current
        #         working directory, although each component
        #         can maintain its own input and output
        #         directories.
        #------------------------------------------------------
        if (self.cfg_prefix == None):
            if (self.USE_GUI_SETTINGS):
                self.cfg_prefix = self.dialog.getString('case_prefix', 'ERROR')
            else:
                print '###########################################'
                print ' ERROR: cfg_prefix has not been set.'
                print '###########################################'
                print ' '
                
        #------------------------------------------------------
        # Note: A run_model() call or a driver's initialize()
        #       call calling initialize_config_vars() will set
        #       CWD to the location of CFG files.
        #------------------------------------------------------
        # Note: If directories and prefixes are not set in
        #       initialize_config_vars(), then they will
        #       default to CWD and cfg_prefix.
        #------------------------------------------------------
##        cfg_extension = self.get_cfg_extension()
##        filename      = self.cfg_prefix + cfg_extension
##        self.cfg_file = os.path.join( os.getcwd(), filename )
      
        if (self.USE_GUI_SETTINGS):
            #------------------------------------------
            # Load input variables that were obtained
            # from a tabbed dialog in the CMT tool.
            #------------------------------------------
            self.read_config_gui()

            #-------------------------------------------
            # Save GUI settings into a new CFG file ?
            #-------------------------------------------
            # Don't need this because all settings can
            # be saved to and loaded from a BLD file.
            #-------------------------------------------
            # self.write_config_file()    AND/OR
            # self.write_gui_info_file()
        else:       
            #---------------------------------------
            # Read input variables from a CFG file
            #---------------------------------------
            # CFG filename was built and saved in
            # the initialize() method before now.
            #---------------------------------------
            ### print '####### CALLING read_config_file()...'
            self.read_config_file()  # (not quite ready yet)
            ### self.read_cfg_file()

        # print 'FINISHED with read_config_gui() ########\n'
        
        #-------------------------------------------------
        # The driver gets to set CWD to its in_directory
        #-------------------------------------------------
        # print 'CHECKING driver mode...'
        if (self.mode == "driver"):
            if (hasattr(self, 'in_directory')):
                if (self.in_directory != None):
                    os.chdir( self.in_directory )

        #-------------------------------------
        # Check the directories and prefixes
        #-------------------------------------
        # Defaults are set in __init__() and
        # some may have just been read in.
        #-------------------------------------
        # print 'CALLING check_directories()...'
        self.check_directories()
                
        #-------------------------------------------------------
        # If component is in "driver" mode, then set its
        # in_directory as current working directory (CWD).
        #
        # There are several reasons for this:
        #
        # (1) Each component's "read_cfg_file()" method
        #     needs to know where to find its CFG file.
        #     We _could_ store the paths to all of a project's
        #     CFG files in the sync_file (created by driver).
        #     However, assuming that all of the CFG files are
        #     in the CWD is simpler and keeps them together.
        #
        # (2) If a user is running an example w/o the GUI,
        #     then it is simplest if the CFG files needed for
        #     the example are stored together with the input
        #     data for the example in a central place on
        #     beach.  It is likely that the user won't have
        #     write permission there, which is why we need a
        #     separate output directory for model output.
        #     We also want to avoid the need for users to
        #     create copies of example CFG files in their own
        #     directories. (But they will have read permission
        #     and can make copies if they want.)
        #
        # (3) Component CFG files represent model input and
        #     may contain input parameters and/or the names
        #     of input files, such as initial-value grids.
        #     If a user is not running an example, then they
        #     will need to create an appropriate set of input
        #     files using whatever tools we give them.
        #
        # (4) Except when running examples (or using someone
        #     else's input files) the directories for input
        #     and output will typically be the same.
        #
        #-------------------------------------------------------
        # The driver also creates a "sync file" to share
        # the site_prefix, case_prefix, etc. with
        # the nondrivers.  This needs to be a name that is
        # not already in use and may be stored in "/tmp".
        #
        # self.mode is set to "driver" or "nondriver" in the
        # initialize() method.
        #-------------------------------------------------------
##        if (self.mode == "driver"):
##            os.chdir( self.in_directory )
##            ## self.write_sync_file()
##        else:
##            #------------------------------------------------
##            # Use information from sync file, such as
##            # site_prefix and case_prefix, to override
##            # those of a component in "nondriver" mode.
##            #------------------------------------------------
##            # The driver passes one argument, the name of
##            # the sync_file to the initialize() method of
##            # each nondriver component.
##            #------------------------------------------------            
##            ## self.read_sync_file()
##            pass
        
        #--------------------------------------------
        # Set any input variables that are computed
        #--------------------------------------------
        # print 'CALLING set_computed_input_vars()...'
        self.set_computed_input_vars()
        
        #-----------------------------------------------------
        # Not all components need this, so don't do it here.
        #-----------------------------------------------------
        # self.read_grid_info()   # (stores rti in self, adds "da")
        
    #   initialize_config_vars()
    #-------------------------------------------------------------------
##    def open_new_sync_file(self):
##
##        #-----------------------------------------------
##        # Make sure this name is not already in use.
##        #
##        # Also need to delete old sync files from /tmp
##        # so they don't accumulate there.
##        #-----------------------------------------------
##        # import os.path
##        # sync_base = '/tmp/project_sync_file_'
##        #-----------------------------------------------        
##        # sync_file = sync_base + '1.txt'
##        # k = 1
##        # while (os.path.exists( sync_file )):
##        #     k += 1
##        #     sync_file = sync_base + str(k) + '.txt'
##        # self.sync_file = sync_file
##        #-----------------------------------------------       
##        # k = 0
##        # while (True):
##        #     k += 1
##        #     sync_file = sync_base + str(k) + '.txt'
##        #     if not(os.path.exits( sync_file )):
##        #         break
##        # self.sync_file = sync_file
##        #-----------------------------------------------        
##        # self.sync_file = '/tmp/project_sync_file.txt'
##        # self.sync_file = '/tmp/driver_sync_file.txt'
##        self.sync_file = 'driver_sync_file.txt'   # (in CWD)
##        
##        file_unit = open( self.sync_file, 'w' )
##        return file_unit
##    
##    #   open_new_sync_file()
##    #-------------------------------------------------------------------
##    def write_sync_file(self, REPORT=False):
##
##        #-------------------------------------------------------
##        # Notes: A "sync file" is a text file with key-value
##        #        pairs that a driver component writes to the
##        #        current working directory (CWD).  It contains
##        #        information such as a site_prefix and a
##        #        case_prefix that are shared by a set of
##        #        coupled components.
##        #
##        #        Drivers call write_sync_file() in their
##        #        initialize() method.  NonDrivers call
##        #        read_sync_file() using a sync_file name that
##        #        the Driver passes to them as the argument to
##        #        their initialize() method.
##        #-------------------------------------------------------
##        file_unit = self.open_new_sync_file()
##
##        if (REPORT):
##            print 'Writing sync_file:', self.sync_file
##            
##        #-------------------------------------------------
##        # Write any information that the Driver wants to
##        # share with the NonDrivers as key-value pairs.
##        #-------------------------------------------------
##        file_unit.write('input_directory:  ' + self.in_directory  + '\n')
##        file_unit.write('output_directory: ' + self.out_directory + '\n')
##        file_unit.write('site_prefix:      ' + self.site_prefix   + '\n')
##        file_unit.write('case_prefix:      ' + self.case_prefix   + '\n')
##        # file_unit.write('some_float:       ' + '3.14159' + '\n')
##
##        file_unit.close()
##            
##    #   write_sync_file()
##    #-------------------------------------------------------------------
##    def read_sync_file(self, REPORT=False):
##
##        if (REPORT):
##            print 'Reading sync_file:', self.sync_file
##            
##        file_unit = open(self.sync_file, 'r')
##
##        #-------------------------
##        # Skip over header lines
##        #-------------------------
##        # cfg.skip_header( file_unit, n_lines=4 )
##
##        while (True):
##            #-------------------------------------
##            # Read a key-value pair as 2 strings
##            #-------------------------------------
##            key, value_str = cfg.read_key_value_pair(file_unit)
##            if (key == ''): break
##
##            if (REPORT):
##                print 'key, value = ' + key + ', ' + value_str
##            
##            #------------------------------------
##            # Save value in self using key name
##            #------------------------------------
##            try:
##                value = numpy.float(value_str)
##                exec('self.' + key + ' = value')
##            except:
##                exec('self.' + key + ' = "' + value_str + '"')
##
##        file_unit.close()
##        
##    #   read_sync_file()  
    #-------------------------------------------------------------------
    def set_computed_input_vars(self):

        pass
    
    #   set_computed_input_vars()
    #-------------------------------------------------------------------
    def initialize_required_components(self, mode="nondriver"):

##        print "MODE     =", mode
##        print "self.CCA =", self.CCA
##        print ' '
        
        if (mode == "driver"):
            if not(self.CCA):
                #-----------------------------------------
                # Create embedded "process object ports"
                #-----------------------------------------
                self.embed_child_components()
                self.add_child_ports()
            if not(self.SILENT):
                print '---------------------------------------'
                ## print 'Initializing CCA ports...'                
                print 'Initializing required components...'
                print '---------------------------------------'
            self.initialize_ports()
        elif (mode == "nondriver"):
            #-----------------------------------------------
            # In "initialize()" method of component's Impl
            # file it must call "get_ports() before
            # calling "initialize()".
            #-----------------------------------------------
            pass
        
    #   initialize_required_components()
    #-------------------------------------------------------------------
    def get_cca_ports(self, d_services):

        #-------------------------------------------------------------
        # Example call in initialize method of "TopoFlow_Impl.py":
        #
        # OK = self.tf.get_cca_ports( self.d_services )
        # if not(OK): return -1
        #-------------------------------------------------------------
        # print 'CALLING get_cca_ports()...'
        
        self.get_cca_port_info()
        port_names   = self.cca_port_names
        short_names  = self.cca_port_short
        port_type    = self.cca_port_type
        project_name = self.cca_project

        exec("import " + project_name)
        
        #--------------------------------------------
        # Does this component have any uses ports ?
        #--------------------------------------------
        if (port_names[0] == ''):
            if (self.DEBUG):
                print '------------------------------------------'
                print 'NOTE: This component has no uses ports.'
                print '------------------------------------------'
            SUCCESS = True
            return SUCCESS
        
        #--------------------------------------------
        # Note:  Use "d_services" added by Bocca to
        #        an Impl file to get a CCA port
        #--------------------------------------------
        SUCCESS = True
        str2 = project_name + "." + port_type + "." + port_type
        str3 = "( port )"
        
        for k in xrange(len(port_names)):
            #-----------------------------------
            # Try to get a "generic" CCA port.
            #-----------------------------------
            name     = port_names[k]
            name_str = '"' + name + '".'
            try:
                port = d_services.getPort( name )
                print "SUCCESS: Got CCA port: " + name_str
                ## print "*** type(port) =", type(port)
            except:
                print '#####################################################'
                print ' ERROR:  Unable to get CCA port named: ' + name_str
                print ' '
                print ' Please select a component from the palette that'
                print ' provides a port named ' + name_str
                print '#####################################################'
                print ' '
                ## message = "FAILURE:  Unable to get CCA port: " + name
                ## print message
                port    = None
                SUCCESS = False

            if not(SUCCESS):
                break

            #------------------------------------------
            # Try to typecast the port to "port_type"
            # and then store it within "self"
            #------------------------------------------
            str1 = "self." + short_names[k]  
            exec( str1 + " = " + str2 + str3 )
 
            exec( "UNABLE = (" + str1 + " == None)" )
            if (UNABLE):
                print 'FAILURE: Unable to cast CCA port: ' + name
                exec( "d_services.releasePort( "  + name + " )")
                SUCCESS = False
           
        return SUCCESS
        
    #   get_cca_ports()
    #-------------------------------------------------------------------
    def release_cca_ports(self, d_services):
   
        #----------------------
        #  Release all ports
        #----------------------
        for name in self.cca_port_names:
            d_services.releasePort( name )
                 
    #   release_cca_ports()    
    #-------------------------------------------------------------------
    def add_child_port(self, port1_str, port2_str, SELF=False):

        #-----------------------------------
        # Example:  port1_str = "cp"
        #           port2_str = "mp"
        #     =>    self.cp.mp = self.mp
        #--------------------------------------
        # Note: If (self == channels (cp)),
        #       then want "self.ep.cp = self"
        #--------------------------------------
        str1 = "self." + port1_str + "."
        str2 = port2_str + " = "
        if not(SELF):
            str3 = "self." + port2_str
        else:
            str3 = "self"
        exec( str1 + str2 + str3)
            
    #   add_child_port()
    #-------------------------------------------------------------------
    def get_port_data(self, var_name, port=None,
                      port_name='CCA', VERBOSE=False):
        
        #---------------------------------------------------------
        # Note: This method is "private", or not part of the
        #       exposed component interface (or port), so its
        #       return type (in Python) can be dynamic.
        #       However, the functions called by this one below
        #       are "port functions" and therefore have static
        #       return types.
        #---------------------------------------------------------
        # Note: To maximize performance, we may want to just get
        #       the type for each var_name once during the
        #       initialize() call and save it internally.
        #---------------------------------------------------------
        if (port is None):
            print "ERROR: get_port_data() needs CCA port argument."
            return numpy.float64(0)
        
        #---------------------------------------------------------
        try:
            if (port.is_scalar(var_name)):
                data = port.get_scalar_double( var_name )
                str2 = 'SCALAR'
            elif (port.is_grid(var_name)):
                data = port.get_grid_double( var_name )
                str2 = 'GRID'
            elif (port.is_vector(var_name)):
                data = port.get_vector_double( var_name )
                str2 = 'VECTOR'
            else:
                print 'ERROR in CSDMS_base.get_port_data().'
                print '    ' + var_name + ' is not SCALAR, VECTOR or GRID.'
                print '    Returning scalar value of 0.'
                print ' '
                data = numpy.float64(0)
                str2 = 'UNKNOWN_TYPE'
        except:
            data = numpy.float64(0)
            str2 = 'SCALAR'
            print 'ERROR in CSDMS_base.get_port_data().'
            pn_str = ' from ' + port_name + ' port.'
            print '   Could not get ' + var_name + pn_str
            print '   Returning scalar value of 0.'
            print ' '
             
        #--------------
        # For testing
        #--------------
        if (VERBOSE):
            str1 = '>>> ' + var_name + ' is a '
            str3 = '. type(' + var_name + ') =', type(data)
            print (str1 + str2 + str3)
             
        #----------------------------------
        # Make sure return type is double
        #----------------------------------
        return numpy.float64(data)
        
    #   get_port_data()
    #-------------------------------------------------------------------
    def get_port_data_OLD(self, var_name, port=None):
        
        #---------------------------------------------------------
        # Note: This method is "private", or not part of the
        #       exposed component interface (or port), so its
        #       return type (in Python) can be dynamic.
        #       However, the functions called by this one below
        #       are "port functions" and therefore have static
        #       return types.
        #---------------------------------------------------------
        # Note: To maximize performance, we may want to just get
        #       the type for each var_name once during the
        #       initialize() call and save it internally.
        #---------------------------------------------------------
        if (port is None):
            print "ERROR: get_port_data() needs CCA port argument."
            return -1
        #---------------------------------------------------------        
##        if (port is None):
##            exec("return self." + var_name)
##        else:
        #---------------------------------------------------------
        if (port.is_scalar(var_name)):
            # print '>>> ' + var_name + ' is a SCALAR.'
            ## return port.get_scalar_double( var_name )
            scalar = port.get_scalar_double( var_name )
            # print '>>> type(' + var_name + ') =', type(scalar)
            return numpy.float64( scalar )
        elif (port.is_grid(var_name)):
            # print '>>> ' + var_name + ' is a GRID.'
            ## return port.get_grid_double( var_name )
            grid = port.get_grid_double( var_name )
            # print '>>> type(' + var_name + ') =', type(grid)
            return numpy.float64( grid )
        else:
##            nv = port.get_size( var_name )
##            if (nv == 1):
##                return numpy.float64( 
            print 'ERROR in CSDMS_base.get_port_data().'
            print '    Variable is not SCALAR or GRID.'
            print '    rank(var) =', port.get_rank( var_name )
            print '    Returning 0.'
            print ' '
            return numpy.float64(0)
            
    #   get_port_data_OLD()
    #-------------------------------------------------------------------
    def set_port_data(self, var_name, var, port=None):
        
        #---------------------------------------------------------
        # Note: This method is "private", or not part of the
        #       exposed component interface (or port), so its
        #       return type (in Python) can be dynamic.
        #       However, the functions called by this one below
        #       are "port functions" and therefore have static
        #       return types.
        #---------------------------------------------------------
        # Note: To maximize performance, we may want to just get
        #       the type for each var_name once during the
        #       initialize() call and save it internally.
        #---------------------------------------------------------
        # Note: (2/8/10) Bug fix.  A variable such as "h_swe"
        #       may be initialized to a scalar (in snow comp.)
        #       so "port.is_scalar()" will return True.  But if
        #       "var" is a grid, then we won't be allowed to
        #       "set it" as one.  Solution is to replace:
        #       if (port.is_scalar(var_name)) to:
        #       if (size(var) == 1):
        #---------------------------------------------------------            
##        if (port is None):
##            exec("return self." + var_name)
##        else:
        #----------------------------------------------
        var = numpy.float64( var )   # (5/20/10)
        if (size(var) == 1):
            port.set_scalar_double(var_name, var)
        else:
            port.set_grid_double(var_name, var)
        #----------------------------------------------
##        if port.is_scalar(var_name):
##            port.set_scalar_double(var_name, var)
##        else:
##            port.set_grid_double(var_name, var)

    #   set_port_data()
    #-------------------------------------------------------------------
    def get_rank(self, var_name):

        exec("rank = numpy.rank(self." + var_name + ")")
        return rank
    
    #   get_rank()
    #-------------------------------------------------------------------
    def get_size(self, var_name):

        #-------------------------------------------------------
        # Notes: This is used by a caller to determine the
        #        number of elements a given variable has.
        #        This information can then be used to size
        #        an array, for example.  See "get_rank()".

        #        In a dynamically-typed language like Python,
        #        the dynamic typing can be used with NumPy to
        #        allow very flexible input types.

        # NB!    Right now, case in var_name must be an exact
        #        match.
        #-------------------------------------------------------
        exec("n = numpy.size(self." + var_name + ")")
        return n
    
    #   get_size()
    #-------------------------------------------------------------------
    def read_grid_info(self):

        #------------------------------------------
        # Read grid info from an RTI file that is
        # in the current working directory.
        #------------------------------------------
        if (self.DEBUG):
            print 'Process component: Reading grid info...'
        self.grid_info_file = (self.in_directory +
                               self.site_prefix + '.rti')

        info = rti_files.read_info( self.grid_info_file )
        
        #----------------------
        # Convenient synonyms
        #-----------------------------------------------------
        # Note that "info" has additional derived attributes
        # such as: n_pixels, bpe, grid_size and SWAP_ENDIAN.
        #-----------------------------------------------------
        self.rti = info
        self.nx  = info.ncols
        self.ny  = info.nrows
        
        #------------------------------------------------
        # Get grid cell areas, "da", which is either a
        # scalar (if same for all grid cells) or a grid
        # with default units of "m^2".
        #------------------------------------------------
        self.da = pixels.get_da( info )
        
##        rti = rti_file.rti_file(self.grid_info_file)
##        OK = rti.read_file()
##        if not(OK):
##            print 'ERROR: Could not read grid info from file:'
##            print self.grid_info_file
##            print ' '
##            return
##        self.rti      = rti
##        self.nx       = rti.ncols
##        self.ny       = rti.nrows
##        self.n_pixels = rti.n_pixels
##        self.da = pixels.get_da(rti)
 
    #   read_grid_info()
    #-------------------------------------------------------------------
    def initialize_basin_vars(self):

        import basins
        self.bp = basins.basins_component()

##        print 'self.cfg_prefix  =', self.cfg_prefix
##        print 'self.site_prefix =', self.site_prefix    ##########
##        print 'self.case_prefix =', self.case_prefix
        
        self.bp.site_prefix = self.site_prefix  ######
        
        self.bp.initialize( cfg_prefix=self.cfg_prefix, 
                            SILENT=not(self.DEBUG) )

        #-------------------------------
        # Store the outlet IDs in self
        #-------------------------------
        outlet_IDs = self.bp.outlet_IDs
        outlet_ID  = outlet_IDs[0]
        self.outlet_IDs = (outlet_IDs / self.nx, outlet_IDs % self.nx)
        self.outlet_ID  = (outlet_ID  / self.nx, outlet_ID  % self.nx)
        
##        self.outlet_IDs = outlet_IDs   # (long-int calendar indices)
##        self.outlet_ID  = outlet_ID

        #--------------------------------------------------
        # Before 5/14/10, get outlet_IDs from the
        # basins port of a Basins component.  Also,
        # every initialize() called "store_outlet_IDs()".
        #--------------------------------------------------
        # outlet_IDs = self.bp.get_vector_long('outlet_IDs')
        
    #   initialize_basin_vars()
    #-------------------------------------------------------------------
##    def store_outlet_IDs(self):
##        
##        outlet_IDs = self.bp.get_vector_long('outlet_IDs')
##        outlet_ID  = outlet_IDs[0]
####        self.outlet_IDs = outlet_IDs   # (long-int calendar indices)
####        self.outlet_ID  = outlet_ID
##        self.outlet_IDs = (outlet_IDs / self.nx, outlet_IDs % self.nx)
##        self.outlet_ID  = (outlet_ID  / self.nx, outlet_ID  % self.nx)
##
##    #   store_outlet_IDs()   
    #-------------------------------------------------------------------
    def initialize_time_vars(self, units='seconds'):

        #------------------
        # Start the clock
        #------------------
        self.start_time = time.time()
        
        #--------------------------------
        # Initialize the time variables
        #--------------------------------
        self.time_units = units.lower()
        self.time_index = int32(0)
        self.time       = float64(0)
        self.DONE       = False
        
        #--------------------------
        # Time conversion factors
        #--------------------------
        self.sec_per_year = float64(365) * 24 * 3600
        self.min_per_year = float64(365) * 24 * 60
        
        #-------------------------------------------
        # For backward compatibility with TopoFlow
        #-------------------------------------------
        self.time_sec = float64(0)
        self.time_min = float64(0)
            
        #--------------------------------------------
        # For print_time_and_value() function below
        #--------------------------------------------
        # Substract 100 seconds so we'll always
        # print values at time zero. (6/29/10)
        #--------------------------------------------
        self.last_print_time = time.time() - 100.0
        
##        self.last_check_time  = time.time()  # (for user interrupt)
##        self.last_plot_time   = float32(0.0)   ### CHECK ###
        
    #   initialize_time_vars()
    #-------------------------------------------------------------------
    def update_time(self):

        #---------------------
        # Increment the time
        #---------------------
        self.time_index += 1
        self.time       += self.dt  # (use same units as dt)
        
        if (self.time_units == 'seconds'):
            self.time_sec = self.time                    # [seconds]
            self.time_min = self.time_sec / float64(60)  # [minutes]
        elif (self.time_units == 'years'):
            #-----------------------------------
            # Used by GC2D and Erode (12/4/09)
            #-----------------------------------
            self.time_sec = self.time * self.sec_per_year  ####
            self.time_min = self.time_sec / float64(60)  # [minutes]
            
    #   update_time()
    #-------------------------------------------------------------------
    def print_time_and_value(self, var, var_name='Q_out',
                             units_name='[m^3/s]',
                             interval=2.0,
                             PRINT_INDEX=False):

        #-----------------------------------------
        # (8/2/10) Print message about interval.
        #-----------------------------------------
        if (self.time_index == 0):
            print 'Will print values every', interval, 'seconds.'
            
        #---------------------------------------------------
        # Note: Print the model time, in minutes, and the
        #       current value of "var", at the specified
        #       real-time "interval" (in seconds).
        #---------------------------------------------------
        # Note: Plotting hydrograph at same interval is
        #       generally too infrequent.
        #---------------------------------------------------
        #  self.Tstr is set in TF by initialize_stop_vars()
        #---------------------------------------------------
        elapsed_time = (time.time() - self.last_print_time)
        if (elapsed_time > interval):
            if (self.time_units == 'seconds'):
                cur_time = self.time_min
                time_units_str = ' [min]'
            else:
                cur_time = self.time
                time_units_str = ' [' + self.time_units + ']' 
            time_str = 'Time = ' + ("%10.2f" % cur_time)
            time_str = time_str + time_units_str
            #-------------------------------------------------
            var_str  = var_name + ' = ' + ("%10.5f" % var)
            var_str  = var_str  + ' ' + units_name          
            #-------------------------------------------------      
            print (time_str + ',  ' + var_str)
            #-----------------------------------------------------
            if (PRINT_INDEX):
                print 'n =', self.time_index, 'of', self.n_steps
            #-----------------------------------------------------                
            self.last_print_time = time.time()

    #   print_time_and_value()    
    #-------------------------------------------------------------------
    def print_run_time(self, proc_name='component',
                       sec_digits=4, seconds=None,
                       SILENT=None):
                       ### SUB_PROCESS=False) 

        #------------------------------------------------------
        # If "seconds" argument is only provided for testing.
        # You can provide this value to make sure that the
        # minuts, hours, days, etc. are computed correctly.
        #------------------------------------------------------
        if (seconds == None):    
            finish  = time.time()
            seconds = (finish - self.start_time)

        #----------------------------------
        # Compute minutes, hours and days
        #----------------------------------
        dec_part  = (seconds % float32(1.0))     #(Save decimal part)
        days      = int32(seconds) / int32(86400)
        secs_left = int32(seconds) % int32(86400)
        hours     = (secs_left / int32(3600))
        secs_left = (secs_left % int32(3600))
        minutes   = (secs_left / int32(60))
        seconds   = (secs_left % int32(60))
        #-----------------------------------------
        #hours     = long(seconds)  /  3600L
        #secs_left = long(seconds) mod 3600L
        #minutes   = (secs_left  /  60L)
        #seconds   = (secs_left mod 60L)
        
        #----------------------------
        # Construct the time string
        #----------------------------
        time_string = ''
        #--------------------------------------------------------
        if (days > 0):    
            if (days > 1):    
                e0 = ' days, '
            else:    
                e0 = ' day, '
            time_string += str(days) + e0
        #--------------------------------------------------------
        if (hours > 0):    
            if (hours > 1):    
                e1 = ' hours, '
            else:    
                e1 = ' hour, '
            time_string += str(hours) + e1
        #--------------------------------------------------------
        if (minutes > 0):    
            if (minutes > 1):    
                e2 = ' minutes, '
            else:    
                e2 = ' minute, '
            time_string += str(minutes) + e2
        
        #-----------------------------------------
        # Default is 4 digits after the decimal.
        #-----------------------------------------
        dec_pastr = ('.' + str(dec_part)[2:2+sec_digits])
        time_string += str(seconds) + dec_pastr + ' seconds.'
        
        if not(SILENT):
            print ('Run time for ' + proc_name + ' = ')
            print  time_string
            print ' '
            
##            if (SUB_PROCESS):    
##                PART1 = '>> '
##            else:    
##                PART1 = ''
##            print (PART1 + 'Run time for ' + procname + ' = ')
##            print (PART1 + time_string)
##            print ' '
     
    #   print_run_time()
    #-------------------------------------------------------------------    
    def print_final_report(self, comp_name='Processs component'):

        if (self.mode == 'nondriver'):
            print comp_name + ': Finished.'
            return
        
        #-------------------
        # Print the report
        #-------------------
        hline = ''.ljust(60, '-')
        print hline
        print tf_utils.TF_Version()      # (specific to TopoFlow)
        print time.asctime()
        print ' '
        print 'Input directory:      ' + self.in_directory
        print 'Output directory:     ' + self.out_directory
        print 'Site prefix:          ' + self.site_prefix
        print 'Case prefix:          ' + self.case_prefix
        print ' '

        #----------------------------
        # Construct run time string
        #----------------------------
        run_time = (time.time() - self.start_time)
        if (run_time > 60):
            run_time = run_time / numpy.float64(60)
            rt_units = ' [min]'
        else:
            rt_units = ' [sec]'
        run_time_str = str(run_time) + rt_units
        
        print 'Simulated time:      ' + str(self.time_min) + ' [min]'
        print 'Program run time:    ' + run_time_str
        print ' '
        print 'Number of timesteps: ' + str(self.time_index)
        print 'Process timestep:    ' + str(self.dt) + ' [s]'  #######
        print 'Number of columns:   ' + str(self.nx)
        print 'Number of rows:      ' + str(self.ny)
        print ' '
        finish_str = ': Finished. (' + self.case_prefix + ')'
        print comp_name + finish_str
        print ' '
        
    #   print_final_report()
    #-------------------------------------------------------------------    
    def print_traceback(self, caller_name='TopoFlow'):

        print '################################################'
        print ' ERROR encountered in ' + caller_name
        print '       Please check your input parameters.'
        print '################################################'
        
        traceback.print_exc()
        
    #   print_traceback()                    
    #-------------------------------------------------------------------

                           
