import QtQuick 2.8
import JASP.Controls 1.0
import JASP.Theme 1.0
import JASP.Widgets 1.0

Form
{
  VariablesForm
  {
    height: 300
    AvailableVariablesList { name: "allVariablesList" }		
    AssignedVariablesList { name: "variables"; title: qsTr("Variables")}
  }
  
  Section
  {
    title: qsTr("Single-Test Reliability")
    Group
    {
      title: qsTr("Scale Statistics")
      CheckBox 
      {    
        name:   "mcDonaldScalef"	
        label:  qsTr("McDonald's ω")         
        id:     mcdonaldf
        
        CheckBox
        {
          name:     "fitMeasures"	
          label:    qsTr("Single Factor Model Fit")         
          enabled:  mcdonaldf.checked
        }
      }
      CheckBox 
      {     
        name: "alphaScalef";				
        label: qsTr("Cronbach's α");         
        id: cronbachf   		      
      }
      CheckBox 
      {     
        name: "guttman2Scalef";			
        label: qsTr("Guttman's λ2");         
        id: guttman2f       	          
      }
//      CheckBox 
  //    {     
//        name: "guttman6Scalef";			
//        label: qsTr("Guttman's λ6");         
//        id: guttman6f       	          
//      }
      CheckBox 
      {     
        name: "glbScalef";				  
        label: qsTr("Greatest lower bound"); 
        id: glbf      	              
      }
      CheckBox { name: "averageInterItemCor";				label: qsTr("Average interitem correlation")	}

      CIField 
      {      
        name: "ConfidenceIntervalValue";   
        label: qsTr("Confidence interval")    
      }
      CheckBox { name: "meanScale";						label: qsTr("Mean")							}
      CheckBox { name: "sdScale";							label: qsTr("Standard deviation")			}
      
    }
    
    Group
    {
      title: qsTr("Individual Item Statistics")
      CheckBox 
      { 
        name: "mcDonaldItemf";				
        label: qsTr("McDonald's ω  (if item dropped)");	        
        enabled: mcdonaldf.checked 
      }
      CheckBox 
      { 
        name: "alphaItemf";				
        label: qsTr("Cronbach's α (if item dropped)");	     
        enabled: cronbachf.checked 
      }
      CheckBox 
      { 
        name: "guttman2Itemf";			
        label: qsTr("Guttman's λ2 (if item dropped)");	       
        enabled: guttman2f.checked  
      }
      CheckBox 
      { 
        name: "glbItemf";     		
        label: qsTr("Greatest lower bound (if item dropped)");	
        enabled: glbf.checked     
      }    
      CheckBox { name: "itemRestCor";						label: qsTr("Item-rest correlation")				}
      CheckBox { name: "meanItem";						label: qsTr("Mean")								}
      CheckBox { name: "sdItem";							label: qsTr("Standard deviation")				}
      
      
    }
    
  }
  Section
  {
    title: qsTr("Reverse-Scaled Items")
    
    VariablesForm
    {
      height: 150
      AvailableVariablesList { name: "normalScaledItems";	 title: qsTr("Normal-Scaled Items"); source: "variables" }
      AssignedVariablesList {  name: "reverseScaledItems"; title: qsTr("Reverse-Scaled Items") }
    }
  }
  
  Section
  {
    title: qsTr("Advanced Options")
    IntegerField
    {
      name: "noSamplesf"
      label: qsTr("No. of bootstrap samples")
      defaultValue: 500
      fieldWidth: 50
      min: 100
      max: 1e7
    }
  RadioButtonGroup {
        title: qsTr("Missing Values")
        name: "missingValues"
        RadioButton { value: "excludeCasesListwise"; label: qsTr("Exclude cases listwise"); checked: true	}
        RadioButton { value: "excludeCasesPairwise"; label: qsTr("Exclude cases pairwise")					}
        }
  }
}
