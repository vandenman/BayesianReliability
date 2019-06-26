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
	
	Group
	{
		title: qsTr("Scale Statistics")
		CheckBox {     name: "mcDonaldScale";			label: qsTr("McDonald's ω");         id: mcdonald;   }
    CheckBox {     name: "alphaScale";				label: qsTr("Cronbach's α");         id: cronbach   		      }
    CheckBox {     name: "guttman2Scale";			label: qsTr("Guttman's λ2");         id: guttman       	          }
		CheckBox {     name: "glbScale";				label: qsTr("Greatest lower bound"); id: glb      	              }

	}
       
	Group
	{
		title: qsTr("Individual Item Statistics")
		CheckBox { name: "mcDonaldItem";				label: qsTr("McDonald's ω  (if item dropped)");	        enabled: mcdonald.checked }
		CheckBox { name: "alphaItem";					label: qsTr("Cronbach's α (if item dropped)");	        enabled: cronbach.checked }
		CheckBox { name: "guttman2Item";				label: qsTr("Guttman's λ2 (if item dropped)");	        enabled: guttman.checked  }
    CheckBox { name: "glbItem";     				label: qsTr("Greatest lower bound (if item dropped)");	enabled: glb.checked      }

        
	}

    Group
    {
        CheckBox {
                       name: "plotPosterior";           label: qsTr("Plot Posteriors");
            CheckBox { name: "dispPrior";               label: qsTr("Display Priors")}         
            CheckBox { name: "fixXRange";               label: qsTr("Fix range to 0-1")}
            CheckBox { name: "cutoff";               label: qsTr("Display cutoffs at"); 
            columns: 2;
              DoubleField
              {
              name: "cutoffValue1"
              label: qsTr("")
              defaultValue: 0.70
              min: 0
              max: 1
              decimals: 2
              fieldWidth: 40
              }               
              DoubleField
              {
              name: "cutoffValue2"
              label: qsTr("")
              defaultValue: 0.80
              min: 0
              max: 1
              decimals: 2
              fieldWidth: 40
              }
            }
        }
        CIField {      name: "CredibleIntervalValue";   label: qsTr("Credible interval")    }
    }
    
    Group
    {
        CheckBox { name: "probTable"; label: qsTr("Probability for Reliability Statistic > x"); 

        DoubleField
        {
            name: "probTableValue"
            label: qsTr("x = ")
            defaultValue: 0.70
            min: 0
            max: 1
            decimals: 2
            fieldWidth: 40

        }
        CheckBox { name: "shadePlots";              label: qsTr("Shade region in plots")     }              
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
            name: "noSamples"
            label: qsTr("No. samples")
            defaultValue: 500
            fieldWidth: 50
            min: 100
            max: 1e7
        }
    }
}
