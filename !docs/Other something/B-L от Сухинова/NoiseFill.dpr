PROGRAM NoiseFill;

{$APPTYPE CONSOLE}

USES
  SysUtils, Graphics;

TYPE
	RGB=Record
		b,g,r: Byte;
	End;
	TPixels=Array[0..0] of RGB;
	PPixels=^TPixels;

	TPoint=Record
		x,y: Integer;
		value: Double;
	End;

	Tlayer=Record
		Frequency: Double;
		Amplitude: Double;
		K: Double;
		Points: Array of TPoint;
	End;

	TLayers=Array of TLayer;

	TField=Record
		Nx,Ny: Integer;
		Min,Max: Single;
		Data: Array of Array of Single;
	End;

Procedure GenPoints(var Layers: TLayers; Nx,Ny: Integer);
Var
	ln: Integer;
	pn: Integer;
	numPoints: Integer;
Begin
	For ln:=0 to High(layers) do Begin
		numPoints:=Round(Nx*Ny*sqr(Layers[ln].Frequency));
		Layers[ln].K:=Nx*Ny/NumPoints;
		SetLength(Layers[ln].Points,numPoints);
		For pn:=0 to High(Layers[ln].Points) do Begin
			Layers[ln].Points[pn].x:=random(Nx);
			Layers[ln].Points[pn].y:=random(Ny);
			Layers[ln].Points[pn].value:=random*Layers[ln].Amplitude;
		End;
	End;
End;

Procedure GenField(const Layers: TLayers; Nx,Ny: Integer; Min,Max: Single; var Field: TField);
Var
	vMin,v,vMax: Double;
	k, ks: Double;
	dx,dy: Double;

	pn,ln,x,y: Integer;
Begin
	If (Nx=0) or (Ny=0) or (Length(Layers)=0) Then Exit;

	SetLength(Field.Data,Nx,Ny);
	Field.Nx:=Nx; Field.Ny:=Ny;

	For y:=0 to Ny-1 do Begin
		For x:=0 to Nx-1 do Begin
			Field.Data[x,y]:=0.0;
		End;
	End;

	For ln:=0 to High(Layers) do Begin
		WriteLn('Processing layer ',ln,'... ');
		For y:=0 to Ny-1 do Begin
			For x:=0 to Nx-1 do Begin
				v:=0.0; ks:=0.0;
				For pn:=0 to High(Layers[ln].Points) do Begin
					dx:=x-Layers[ln].Points[pn].x;
					dy:=y-Layers[ln].Points[pn].y;
					dx:=dx-Round(dx/Nx)*Nx;
					dy:=dy-Round(dy/Ny)*Ny;

					k:=1.0/(sqr(dx)+sqr(dy)+Layers[ln].K);
					v:=v+Layers[ln].Points[pn].value*k;
					ks:=ks+k;
				End;
				Field.Data[x,y]:=Field.Data[x,y]+v/ks;
			End;
		End;
	End;

	vMin:=Field.Data[0,0]; vMax:=vMin;
	For y:=0 to Ny-1 do Begin
		For x:=0 to Nx-1 do Begin
			v:=Field.Data[x,y];
			if v<vMin then vMin:=v;
			if v>vMax then vMax:=v;
		End;
	End;

	For y:=0 to Ny-1 do Begin
		For x:=0 to Nx-1 do Begin
			Field.Data[x,y]:=(Field.Data[x,y]-vMin)/(vMax-vMin)*(Max-Min)+Min;
		End;
	End;
	Field.Min:=Min; Field.Max:=Max;
End;

Procedure ReduceField(var Field: TField);
Var
 x,y: Integer;
 v1,v2,v3,v4: Double;
Begin
  For y:=0 to (Field.Ny div 2)-1 do Begin
    For x:=0 to (Field.Nx div 2)-1 do Begin
      v1:=Field.Data[x*2,y*2  ]; v2:=Field.Data[x*2+1,y*2  ];
      v3:=Field.Data[x*2,y*2+1]; v4:=Field.Data[x*2+1,y*2+1];

      Field.Data[x,y]:=( ln(1/(exp(-v1)+exp(-v3))+1/(exp(-v2)+exp(-v4)))
                       + ln(1/(exp(-v1)+exp(-v2))+1/(exp(-v3)+exp(-v4))) )/2;
    End;
  End;

  Field.Nx:=Field.Nx div 2; Field.Ny:=Field.Ny div 2;
  SetLength(Field.Data, Field.Nx, Field.Ny);
End;

Function DrawField(const Field: TField): TBitmap;
Var
	C: RGB;
	SL: PPixels;

	x,y: Integer;
	v: Single;
Begin
	result:=TBitmap.Create;
	result.PixelFormat:=pf24bit;
	result.Width:=Field.Nx; result.Height:=Field.Ny;

	For y:=0 to Field.Ny-1 do Begin
		SL:=result.ScanLine[y];
		For x:=0 to Field.Nx-1 do Begin
			v:=Field.Data[x,y];
			C.r:=Round((v-Field.Min)/(Field.Max-Field.Min)*255);
			C.g:=C.r;
			C.b:=C.r;
			SL^[x]:=C;
		End;
	End;
End;

Procedure WriteField(const Field: TField; FileName: String);
Var
	x,y: Integer;

	F: Text;
Begin
	Assign(F, FileName);
	ReWrite(F);

	WriteLn(F, Field.Nx,' ',Field.Ny);
	For y:=0 to Field.Ny-1 do Begin
		For x:=0 to Field.Nx-1 do Begin
			Write(F, Field.Data[x,y],' ');
		End;
		WriteLn(F);
	End;

	Close(F);
End;

Var
	Nx,Ny: Integer;
	Nl: Integer;
	Min,Max: Single;

	Layers: TLayers;

	n: Integer;

	Field: TField;
	BM: TBitmap;
begin
  Randomize;
  Write('Enter NX: '); ReadLn(Nx);
  Write('Enter NY: '); ReadLn(Ny);
  Write('Enter NumLayers: '); ReadLn(Nl);
  Write('Enter minimum value: '); ReadLn(Min);
  Write('Enter maximum value: '); ReadLn(Max);

  SetLength(Layers,Nl);

  For n:=0 to Nl-1 do Begin
    WriteLn('Layer ',n);
    Write('Enter 0<Frequency<1: '); ReadLn(Layers[n].Frequency);
    Write('Enter Amplitude: '); ReadLn(Layers[n].Amplitude);
  End;

  WriteLn;
  Write('Generating points... ');
  GenPoints(Layers,Nx,Ny);
  WriteLn('Done.');

  WriteLn;
  WriteLn('Generating field... ');
  GenField(Layers,Nx,Ny, Min,Max, Field);
  WriteLn('Done.');

  WriteLn;
  Write('Saving result to files... ');

  WriteField(Field, 'Res.dat');
  BM:=DrawField(Field);
  BM.SaveToFile('Res.bmp');
  BM.Free;

  ReduceField(Field);

  WriteField(Field, 'Res2.dat');
  BM:=DrawField(Field);
  BM.SaveToFile('Res2.bmp');
  BM.Free;

  WriteLn('Done.');
  ReadLn;
end.
