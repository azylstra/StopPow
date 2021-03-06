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

import cStopPow.StopPow;
import java.util.concurrent.ThreadPoolExecutor;
import javax.swing.SwingUtilities;

/**
 * Generate a JPanel containing UI elements allowing the user
 * to do calculations for stopping power
 * @author Alex Zylstra
 * @class CalculatorPaneldEdx
 * @date 2013/05/18
 */
public class CalculatorPaneldEdx extends CalculatorPanel
{
    /** The calculation result */
    float calcResult;
    
    /**
     * Create a new CalculatorPaneldEdx.
     * @param key The String key for the model that this panel is using. 
     * Must match the key for model in the ModelManager. Can be null
     * when instantiating a new panel.
     * @param manager The ModelManager in use by this application
     * @param calcPool A ThreadPool for use when doing background calculations
     */
    public CalculatorPaneldEdx(String key, ModelManager manager, ThreadPoolExecutor calcPool)
    {
        super(key,manager,calcPool);
        
        initComponents();
        
        // default starting value:
        calcResult = Float.NaN;
    }
    
    /**
     * Update the underlying stopping power calculation. Does not launch in
     * a new thread, use the Runnable interface for that:
     * {@code exec.execute(myCalculatorPanelShift)}
     * where exec is a ThreadPool.
     */
    @Override
    public void updateCalculation() 
    {
        // sanity check:
        if( key == null || !manager.containsKey(key) )
        {
            // if model is not set, or no energy is entered,
            // clear text:
            calcResult = Float.NaN;
            return;
        }
        
        // try to gray out input 
        try {
            energyTextField.setEditable(false);
        }
        catch(Exception e) {}
        
        // need to lock the panel and the model:
        panelLock.lock();
        manager.lockModel(key);
        try
        {   
            // get the model
            StopPow model = manager.get_model(key);
            
            // set the model's mode based on drop-down menu:
            int modeIndex = modeComboBox.getSelectedIndex();
            if( modeComboBox.getItemAt(modeIndex) == "MeV/(mg/cm2)")
                model.set_mode( StopPow.getMODE_RHOR() );
            else if( modeComboBox.getItemAt(modeIndex) == "MeV/um")
                model.set_mode( StopPow.getMODE_LENGTH() );
            
            // get the energy:
            float energy;
            try {
                energy = Float.parseFloat( energyTextField.getText() );
            }
            catch(java.lang.NumberFormatException e) {
                energy = Float.NaN;
            }

            // calculate the stopping power
            // wrap in try/catch in case the StopPow calculation throws
            try {
                calcResult = model.dEdx(energy);
            }
            catch (java.lang.IllegalArgumentException e) {
                calcResult = Float.NaN;
            }
            catch (Exception e) {
                calcResult = Float.NaN;
            }
        }
        finally
        {
            panelLock.unlock();
            manager.unlockModel(key);
        }
        
        // update the GUI:
        // use SwingUtilities to let the EDT do this later:
        SwingUtilities.invokeLater(new CalculatorPanelGUIRunnable(this));
    }
    
    /**
     * Clear the displayed result.
     */
    @Override
    public void clearResults()
    {
        dEdxTextField.setText("");
    }
    
    /**
     * Update the displayed dynamic elements to reflect the latest
     * calculation result and mode selection.
     */
    @Override
    public void updateGUI()
    {
        // make sure the input is editable:
        energyTextField.setEditable(true);
        
        // set the displayed calculation result:
        if( Float.isNaN(calcResult) )
            dEdxTextField.setText( "" );
        else
            dEdxTextField.setText( String.format(outputFormat, calcResult) );
    }
    
    /**
     * This method is called from within the constructor to initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is always
     * regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jLabel1 = new javax.swing.JLabel();
        energyTextField = new javax.swing.JTextField();
        jLabel2 = new javax.swing.JLabel();
        jLabel3 = new javax.swing.JLabel();
        modeComboBox = new javax.swing.JComboBox();
        jLabel4 = new javax.swing.JLabel();
        dEdxTextField = new javax.swing.JTextField();

        jLabel1.setText("Energy: ");

        energyTextField.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                energyTextFieldActionPerformed(evt);
            }
        });

        jLabel2.setText("MeV");

        jLabel3.setText("Mode: ");

        modeComboBox.setModel(new javax.swing.DefaultComboBoxModel(new String[] { "MeV/um", "MeV/(mg/cm2)" }));
        modeComboBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                modeComboBoxActionPerformed(evt);
            }
        });

        jLabel4.setText("dE/dx: ");

        dEdxTextField.setEditable(false);
        dEdxTextField.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseReleased(java.awt.event.MouseEvent evt) {
                dEdxTextFieldMouseReleased(evt);
            }
        });
        dEdxTextField.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                dEdxTextFieldActionPerformed(evt);
            }
        });

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(this);
        this.setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jLabel4)
                        .addGap(20, 20, 20)
                        .addComponent(dEdxTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 90, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jLabel1)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(energyTextField, javax.swing.GroupLayout.PREFERRED_SIZE, 90, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jLabel2))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jLabel3)
                        .addGap(18, 18, 18)
                        .addComponent(modeComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel1)
                    .addComponent(energyTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel2))
                .addGap(20, 20, 20)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(modeComboBox, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                    .addComponent(jLabel3))
                .addGap(20, 20, 20)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(jLabel4)
                    .addComponent(dEdxTextField, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap())
        );
    }// </editor-fold>//GEN-END:initComponents

    private void energyTextFieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_energyTextFieldActionPerformed
        exec.execute(this);
    }//GEN-LAST:event_energyTextFieldActionPerformed

    private void modeComboBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_modeComboBoxActionPerformed
        exec.execute(this);
    }//GEN-LAST:event_modeComboBoxActionPerformed

    private void dEdxTextFieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_dEdxTextFieldActionPerformed
        // 
    }//GEN-LAST:event_dEdxTextFieldActionPerformed

    private void dEdxTextFieldMouseReleased(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_dEdxTextFieldMouseReleased
        // check for right mouse button:
        if( evt.getButton() == java.awt.event.MouseEvent.BUTTON3 )
        {
            // calculate coordinates:
            int x = evt.getX() + dEdxTextField.getX();
            int y = evt.getY() + dEdxTextField.getY();
            rightClickMenu.show(this, x, y);
        }
    }//GEN-LAST:event_dEdxTextFieldMouseReleased

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JTextField dEdxTextField;
    private javax.swing.JTextField energyTextField;
    private javax.swing.JLabel jLabel1;
    private javax.swing.JLabel jLabel2;
    private javax.swing.JLabel jLabel3;
    private javax.swing.JLabel jLabel4;
    private javax.swing.JComboBox modeComboBox;
    // End of variables declaration//GEN-END:variables
}
