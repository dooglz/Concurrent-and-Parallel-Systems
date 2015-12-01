//*********************************************************
//
// Copyright (c) Microsoft. All rights reserved.
// This code is licensed under the MIT License (MIT).
// THIS CODE IS PROVIDED *AS IS* WITHOUT WARRANTY OF
// ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING ANY
// IMPLIED WARRANTIES OF FITNESS FOR A PARTICULAR
// PURPOSE, MERCHANTABILITY, OR NON-INFRINGEMENT.
//
//*********************************************************

#include "stdafx.h"
#include "D3D12nBodyGravity.h"
#include <iostream>
#include <string> 
#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable :4996)

unsigned int ParticleCount = 4096;

_Use_decl_annotations_
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE, LPSTR, int nCmdShow)
{
  AllocConsole();
  AttachConsole(GetCurrentProcessId());
  freopen("CONIN$", "r", stdin);
  freopen("CONOUT$", "w", stdout);
  freopen("CONOUT$", "w", stderr);

  std::cout << __argc << std::endl;
  if (__argc > 1) {
    ParticleCount = (unsigned int)std::stoi(__argv[1]);
  }
  else {
    ParticleCount = 8000;
  }
  std::cout << "Particle count is: " << ParticleCount << std::endl;

  D3D12nBodyGravity sample(1280, 1280, L"D3D12 n-Body Gravity Simulation");
  return Win32Application::Run(&sample, hInstance, nCmdShow);
}
