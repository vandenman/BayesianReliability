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
      title: qsTr("Bayesian Single-Test Reliability")
    	Group
    	{
    		title: qsTr("Scale Statistics")
    		CheckBox 
    		{    
    		  name:   "mcDonaldScale"	
    		  label:  qsTr("McDonald's ω")         
    		  id:     mcdonald
    
    		  CheckBox 
    		  {    
    		  name:   "dispPPC"	
    		  label:  qsTr("Display posterior predictive check"); enabled: mcdonald.checked        
    	  	}
    
    		}
        CheckBox 
        {     
          name: "alphaScale";				
          label: qsTr("Cronbach's α");         
          id: cronbach   		      
          }
        CheckBox 
        {     
          name: "guttman2Scale";			
          label: qsTr("Guttman's λ2");        
          id: guttman      
          }
    		CheckBox 
    		{     
    		  name: "glbScale";				  
    		  label: qsTr("Greatest lower bound"); 
    		  id: glb      	              
    		  }
    		CIField 
    		{      
    		name: "CredibleIntervalValue";   
    		label: qsTr("Credible interval")    
    		}
    
    	}
           
    	Group
    	{
    		title: qsTr("Individual Item Statistics")
    		CheckBox 
    		{ 
    		  name: "mcDonaldItem";				
    		  label: qsTr("McDonald's ω  (if item dropped)");	        
    		  enabled: mcdonald.checked 
    		  }
    		CheckBox 
    		{ 
    		  name: "alphaItem";					
    		  label: qsTr("Cronbach's α (if item dropped)");	        
    		  enabled: cronbach.checked 
    		  }
    		CheckBox 
    		{ 
    		  name: "guttman2Item";				
    		  label: qsTr("Guttman's λ2 (if item dropped)");	        
    		  enabled: guttman.checked  
    		  }
        CheckBox 
        { 
          name: "glbItem";     				
          label: qsTr("Greatest lower bound (if item dropped)");	
          enabled: glb.checked     
          }
    	}
    
      Group
      {
          CheckBox 
          {
            name: "plotPosterior";           
            label: qsTr("Plot Posteriors");
            CheckBox 
            {   
              name: "fixXRange";               
              label: qsTr("Fix range to 0-1")
              }
            CheckBox 
            { 
              name: "cutoff";               
              label: qsTr("Display cutoffs at")
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
              CheckBox { 
                name: "dispPrior";               
                label: qsTr("Display Priors")
                }      
            }
        }
        
        Group
        {
            CheckBox 
            { 
              id:                 probTable
              name:               "probTable"; 
              label:              qsTr("Probability for Reliability Statistic  >")
              childrenOnSameRow:  true
            
              DoubleField
              {
                
                  name: "probTableValue"
                  defaultValue: 0.70
                  min: 0
                  max: 1
                  decimals: 2
                  fieldWidth: 40
      
                }
              }
              Item
              {
                width:  shadePlots.width + Theme.subOptionOffset
                height: shadePlots.height
                
                CheckBox 
                { 
                  id:       shadePlots
                  name:     "shadePlots";              
                  label:    qsTr("Shade region in plots"); 
                  enabled:  probTable.checked    
                  x:        Theme.subOptionOffset
                }
              }
          }
      }
      
      
    Section
    {
      title: qsTr("Frequentist Single-Test Reliability")
    	Group
    	{
    		title: qsTr("Scale Statistics")
    		CheckBox 
    		{    
    		  name:   "mcDonaldScalef"	
    		  label:  qsTr("McDonald's ω")         
    		  id:     mcdonaldf
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
          id: guttmanf       	          
          }
    		CheckBox 
    		{     
    		  name: "glbScalef";				  
    		  label: qsTr("Greatest lower bound"); 
    		  id: glbf      	              
    		  }
    		CIField 
    		{      
    		  name: "ConfidenceIntervalValue";   
    		  label: qsTr("Confidence interval")    
    		  }
    
    
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
    		  enabled: guttmanf.checked  
    		  }
        CheckBox 
        { 
          name: "glbItemf";     		
          label: qsTr("Greatest lower bound (if item dropped)");	
          enabled: glbf.checked     
          }
    
            
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
            label: qsTr("No. of posterior samples")
            defaultValue: 500
            fieldWidth: 50
            min: 100
            max: 1e7
        }
        IntegerField
        {
            name: "noSamplesf"
            label: qsTr("No. of bootstrap samples")
            defaultValue: 500
            fieldWidth: 50
            min: 100
            max: 1e7
        }
    }
}
