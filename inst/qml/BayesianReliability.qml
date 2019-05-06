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
		CheckBox {     name: "mcDonaldScale";			label: qsTr("McDonald's ω");         id: mcdonald;   checked: true}
        CheckBox {     name: "alphaScale";				label: qsTr("Cronbach's α");         id: cronbach   		      }
        CheckBox {     name: "gutmann2Scale";			label: qsTr("Gutmann's λ2");         id: gutmann       	          }
		CheckBox {     name: "glbScale";				label: qsTr("Greatest lower bound"); id: glb      	              }
		CheckBox {     name: "meanScale";				label: qsTr("Mean");                 id: mean                     }
		CheckBox {     name: "sdScale";					label: qsTr("Standard deviation");   id: sd                       }
	}
       
	Group
	{
		title: qsTr("Individual Item Statistics")
		CheckBox { name: "mcDonaldItem";				label: qsTr("McDonald's ω  (if item dropped)");	        enabled: mcdonald.checked }
		CheckBox { name: "alphaItem";					label: qsTr("Cronbach's α (if item dropped)");	        enabled: cronbach.checked }
		CheckBox { name: "gutmann2Item";				label: qsTr("Gutmann's λ2 (if item dropped)");	        enabled: gutmann.checked  }
        CheckBox { name: "glbItem";     				label: qsTr("Greatest lower bound (if item dropped)");	enabled: glb.checked      }
		CheckBox { name: "meanItem";					label: qsTr("Mean (if item dropped)");				    enabled: mean.checked     }
		CheckBox { name: "sdItem";						label: qsTr("Standard deviation (if item dropped)");	enabled: sd.checked       }
        
	}

    Group
    {
        CheckBox {
                       name: "plotPosterior";           label: qsTr("Plot Posteriors");
            CheckBox { name: "fixXRange";               label: qsTr("Fix range to 0-1"); checked: true      }
        }
        CIField {      name: "CredibleIntervalValue";   label: qsTr("Credible interval")    }
    }
    
    Group
    {
        CheckBox {
                   name: "probTable";               label: qsTr("Probability Statistics > x");
        CIField  { name: "probTableValue";          label: qsTr("x")                         }
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
