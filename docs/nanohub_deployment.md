---
title: "Deploying Tools to NanoHUB"
keywords: NanoBio, NCN, NanoHUB
topnav: topnav
hide_sidebar: true
summary: This is a quick cheat sheet based on recent experiences and is by no means a substitute to any of NanoHUB documentations. Use this for quick reference but refer to official documentation for detailed guidance. 
---

NanoHUB facilitates multiple options to enable user interactivity with the tools. The following steps are particularly tuned for NSF funded NCN nodes to deploy tools to NanoHUB Cyber Platform. The steps outlined illustrate one particular approach. These may or may not apply to your use case. Contact NanoHUB support for appropriateness of these instructions to your use case. 

## Set up workspace


## Source Code Directory Structure 

## Build Instructions - Makefile

## Rappture vs Jupyter 


## Rapptureizing a tool

### GUI Rendering 

### Wrapper Scripts 

### Input Parameters 

### Executing the applications 

### Output Plots

## Deployment Workflow

The following linear workflow is just to get a quick sense. Some of the steps are iterative. 

Status progress: Registered &rarr; Created &rarr; Updated &rarr; Installed &rarr; Approved &rarr; Published

* Register the tool on NanoHUB.org
    * A tool administrator has to manually approve this registration. 
* Providing Source Code
    * The source code is maintained on nanoFORGE svn repository 
    * If the source code is open source and in a public git repository a link to the code can be provided during tool registration. 
    * If the tool is closed source or not using a public git repository, a compressed file of the code can be uploaded during or after the tool registration. 
    * The source code is manually deployed on hub environment by an administrator
* Describe the tool metadata which includes following information 
    * A title for the tool
    * Short description  
    * Abstract or a detailed description about the tool
    * Credits
    * Citations - Typically is auto-generated from NanoHUB link, can use a custom one like a seminal paper about the tool.
    * Funding Sponsors
    * References
    * Publications 
* Updating Source Code
    * Changes to the source code should be communicated to the administrators from the tool status page. Once an update is requested the tool status is changed back to Updated. 
    * NanoHUB tool administrators manually re-deploy the modified code. 
    * The status will be changed to Installed. 
* Once the code is tested and working, the tool owner can approve the tool. 
* The last step is to publish the tool for registered users to see and use it.
* An obvious subset of these will need to be repeated for tool updates. 

