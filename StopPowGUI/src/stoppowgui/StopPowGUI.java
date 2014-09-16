// StopPow - a charged-particle stopping power library
// Copyright (C) 2014  Massachusetts Institute of Technology / Alex Zylstra

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package stoppowgui;

import stoppowgui.util.LibraryLoader;
import stoppowgui.util.ToggleButtonWindowListener;
import SciTK.*;

import javax.swing.JToggleButton;
import javax.swing.JButton;
import javax.swing.JLabel;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.util.concurrent.ScheduledThreadPoolExecutor;
import java.util.concurrent.ThreadPoolExecutor;
import javax.swing.AbstractButton;
import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.SwingUtilities;


/**
 * Main window for the stopping power calculator
 * @author Alex Zylstra
 */
public class StopPowGUI extends SciTK_App
{
    private ModelManager models;
    private ModelManagerGUI modelsGUI;
    private PlotManagerGUI plotGUI;
    protected ThreadPoolExecutor exec;
    
    /** Constructor */
    public StopPowGUI()
    {
        // use constructor for superclass SciTK_App with specified size:
        super(200,200);
        NAME = "Stopping Power Calculator";
        AUTHOR = "Alex Zylstra";
        DATE = "June 06, 2013";
        
        String logo_path = "/stoppowgui/resources/logo.png";
        about_icon = new ImageIcon(getClass().getResource(logo_path));
        tray_icon = new ImageIcon(getClass().getResource(logo_path));
        
        // set up the threading:
        int threads = Runtime.getRuntime().availableProcessors();
        if( threads > 1 )
            threads--; // on multi core systems, leave one free
        exec = new ScheduledThreadPoolExecutor( threads );
        
        // set up models and model manager:
        models = new ModelManager();
        modelsGUI = new ModelManagerGUI(this,models);
        
        // set up the plot manager:
        plotGUI = new PlotManagerGUI(this,models,exec);
        
        initUI();
        
        setLayout(new BoxLayout(getContentPane(), BoxLayout.Y_AXIS));
	setDefaultCloseOperation(EXIT_ON_CLOSE);
        appUI();
    }
    
    /** Handle opening a file. Does nothing for this application. */
    @Override
    public void openFile() {}
    
    /** Load the initial UI */
    @Override
    public void appUI()
    {
        // Set a bit of text:
        JLabel instruction = new JLabel(" Stopping power: ");
        instruction.setFont(new Font("Serif", Font.BOLD, 16));
        instruction.setAlignmentX(0.5f);
        add(instruction);

        // model_manager_ button
        JToggleButton button_model_manager = new JToggleButton("Configure Models");
        modelsGUI.addWindowListener( new ToggleButtonWindowListener(button_model_manager) );
        button_model_manager.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent event) {
                // find out if the button is selected:
                AbstractButton button = (AbstractButton) event.getSource();
                boolean selected = button.getModel().isSelected();
                
                if( selected )
                {
                    // selected, show the model manager:
                    modelsGUI.setVisible(true);
                    button.setSelected(true);
                    
                } 
                else if( !selected )
                {
                    // deselected, hide it:
                    modelsGUI.setVisible(false);
                    button.setSelected(false);
                }
            }
        });
        button_model_manager.setAlignmentX(0.5f);
        add(button_model_manager);
        
        // Calculators button
        JButton button_calculators = new JButton("Calculators");
        button_calculators.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent ae) {
                // launch a new calculator:
                Calculator c = new Calculator(null,models,exec);
            }
        });
        button_calculators.setAlignmentX(0.5f);
        add(button_calculators);
        
        // plots button
        JToggleButton button_plots = new JToggleButton("Plots");
        plotGUI.addWindowListener( new ToggleButtonWindowListener(button_plots) );
        button_plots.addActionListener(new ActionListener() {
            public void itemStateChanged(ItemEvent event) {
                // foo
            }

            @Override
            public void actionPerformed(ActionEvent ae) {
                // find out if the button is selected:
                AbstractButton button = (AbstractButton) ae.getSource();
                boolean selected = button.getModel().isSelected();
                
                if( selected )
                {
                    // selected, show the model manager:
                    plotGUI.setVisible(true);
                    button.setSelected(true);
                    
                } 
                else if( !selected )
                {
                    // deselected, hide it:
                    plotGUI.setVisible(false);
                    button.setSelected(false);
                }
            }
        });
        button_plots.setAlignmentX(0.5f);
        add(button_plots);


    }
    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // use our LibraryLoader functionality:
        try
        {
            LibraryLoader.load_JAR("/cStopPow/libcStopPow.jnilib");
        }
        catch(java.io.IOException e)
        {
            DialogError ex = new DialogError(null,"Could not load from native library");
            return;
        }
        
        /*
        // First try to load the native library directly from a file:
        try
        {
            System.loadLibrary("cStopPow"); 
        }
        catch( UnsatisfiedLinkError e )
        {
            //SciTK_Dialog d = new Dialog_Error("Could not load native library.");
            //return;
        }
        
        // the native library might also be packaged:
        try
        {
            //StopPowLibLoader s = new StopPowLibLoader();
            //s.loadStopPow();
            LibraryLoader.load_JAR("/cStopPow/libcStopPow.jnilib");
        }
        catch( UnsatisfiedLinkError e )
        {
            SciTK_Dialog d = new Dialog_Error("Could not load native library.");
            return;
        }
        catch( java.io.IOException e )
        {
            SciTK_Dialog d = new Dialog_Error(" IO Error: " + e.getMessage());
        }*/
        
        // Create a new NouveauAnalysis object and let
        // Java invoke it when appropriate:
        SwingUtilities.invokeLater(new Runnable() {
                @Override
                public void run() {
                        StopPowGUI g = new StopPowGUI();
                        g.setVisible(true);
                }
        });
    }
}
