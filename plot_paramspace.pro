PRO go

	data=mrdfits('looped_results_varyfedd.fits',1)

;	tmp=fltarr(210)
;	rlims=[tmp+0.1,tmp+0.2,tmp+0.3,tmp+0.4,tmp+0.5,tmp+0.6,tmp+0.7,tmp+0.8,tmp+0.9,tmp+1.]
;	rlim=0.
;	WHILE (n_elements(rlim) LT 21000) DO rlim=[rlim,rlims]  
;	rlim=rlim[1:n_elements(rlim)-1] 
;	data.ratio_lim=rlim

;	tmp=fltarr(2100)
;	fedds=[tmp+0.1,tmp+0.2,tmp+0.3,tmp+0.4,tmp+0.5,tmp+0.6,tmp+0.7,tmp+0.8,tmp+0.9,tmp+1.]
;	fedd=0.
;	WHILE (n_elements(fedd) LT 21000) DO fedd=[fedd,rlims]  
;	fedd=fedd[1:n_elements(fedd)-1] 
;	data.mean_fedd=fedd

	Lcuts=dindgen(30)/10.+44.5
  	ratio_lims=dindgen(10)/10.+0.1
  	obsc_life_fracs=dindgen(7)/10.+0.2
  	mean_fedds=dindgen(10.)/10.+0.1


	lcut=46.8
	ratio_lim=0.6
	obsc_life_frac=0.5
	mean_fedd=0.1
;	lcut=45.9
;	ratio_lim=0.7
;	obsc_life_frac=0.6
;	mean_fedd=0.1

	varlcuts=WHERE(data.ratio_lim EQ ratio_lim AND $
				   data.obsc_life_frac EQ obsc_life_frac AND $
				   data.mean_fedd EQ mean_fedd AND $
				   data.mhalo_obsc_med NE -99 AND $
				   data.mhalo_unob_med NE -99)

	varratio=WHERE(data.lcut EQ lcut AND $
				   data.obsc_life_frac EQ obsc_life_frac AND $
				   data.mean_fedd EQ mean_fedd AND $
				   data.mhalo_obsc_med NE -99 AND $
				   data.mhalo_unob_med NE -99)

	varobsclife=WHERE(data.lcut EQ lcut AND $
				  	  data.ratio_lim EQ ratio_lim AND $
				   	  data.mean_fedd EQ mean_fedd AND $
				   	  data.mhalo_obsc_med NE -99 AND $
				      data.mhalo_unob_med NE -99)

	varfedd=WHERE(data.lcut EQ lcut AND $
				  data.ratio_lim EQ ratio_lim AND $
			 	  data.obsc_life_frac EQ obsc_life_frac AND $
				  data.mhalo_obsc_med NE -99 AND $
				  data.mhalo_unob_med NE -99)

	PS_start,filename='paramspace.eps',xsize=18,ysize=14,/encapsul
 	circsym,/fill
	xtit=textoidl('Log(L_{cut} [ergs/s])')
	ytit1=textoidl('Log(M_{h} [h^{-1} M'+sunsymbol()+'])')
	ytit2=textoidl('b_q')
	ytit3=textoidl('f_{obsc}')
	plot,[0],[0],psym=0,$
        ytit=ytit1,yra=[11.8,13.5],ysty=1,$
        xra=[45.5,47.2],xsty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.1,0.67,0.32,0.95],xtickname=[' ',' ',' ',' ',' ',' ',' ']				   
	oplot,data[varlcuts].lcut,data[varlcuts].mhalo_obsc_med,linestyle=0,color=cgcolor('red'),thick=8
	oplot,data[varlcuts].lcut,data[varlcuts].mhalo_unob_med,linestyle=0,color=cgcolor('blue'),thick=8
	oplot,data[varlcuts].lcut,data[varlcuts].mhalo_obsc_mean,linestyle=2,color=cgcolor('red'),thick=8
	oplot,data[varlcuts].lcut,data[varlcuts].mhalo_unob_mean,linestyle=2,color=cgcolor('blue'),thick=8
	legend,['DiPompeo 16','Mendez 15','Shen 09'],psym=[8,7,2],$
		   color=[cgcolor('black'),cgcolor('purple'),cgcolor('magenta')],$
	       box=0,/bottom,/right, charsize=1.8,charthick=1.8
	
	oplot,[46.84],[13.03],psym=2,symsize=1.5,thick=2,color=cgcolor('magenta')
	oploterror,[46.84],[13.03],0.11,psym=3,thick=5,color=cgcolor('magenta')
	oplot,[46.39],[12.57],psym=2,symsize=1.5,thick=2,color=cgcolor('magenta')
	oploterror,[46.39],[12.57],0.11,psym=3,thick=5,color=cgcolor('magenta')
	
	oplot,[44.5+ALOG10(12)],[11.58],psym=7,symsize=1.5,thick=2,color=cgcolor('purple')
	oploterror,[44.5+ALOG10(12)],[11.58],2.01,psym=3,thick=5,color=cgcolor('purple'),/lobar
	oploterror,[44.5+ALOG10(12)],[11.58],0.72,psym=3,thick=5,color=cgcolor('purple'),/hibar

	circsym	
	oplot,[46.2],[12.94],psym=8,symsize=1.5,thick=2,color=cgcolor('red')
	oploterror,[46.2],[12.94],0.18,psym=3,thick=5,color=cgcolor('red'),/lobar
	oploterror,[46.2],[12.94],0.15,psym=3,thick=5,color=cgcolor('red'),/hibar
	circsym,/fill
	oplot,[46.2],[12.56],psym=8,symsize=1.5,thick=2,color=cgcolor('blue')
	oploterror,[46.2],[12.56],0.21,psym=3,thick=5,color=cgcolor('blue'),/lobar
	oploterror,[46.2],[12.56],0.17,psym=3,thick=5,color=cgcolor('blue'),/hibar
	oplot,[46.2],[13.42],psym=8,symsize=1.5,thick=2,color=cgcolor('red')
	oploterror,[46.2],[13.42],0.39,psym=3,thick=5,color=cgcolor('red')
	
	
	plot,[0],[0],psym=0,$
        ytit=ytit2,yra=[1.1,3.9],ysty=1,$
        xra=[45.5,47.2],xsty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.1,0.39,0.32,0.67],/noerase,xtickname=[' ',' ',' ',' ',' ',' ',' ']		
	oplot,data[varlcuts].lcut,data[varlcuts].b_obsc,linestyle=0,color=cgcolor('red'),thick=8
	oplot,data[varlcuts].lcut,data[varlcuts].b_unob,linestyle=0,color=cgcolor('blue'),thick=8
	
	plot,[0],[0],psym=0,$
        xtit=xtit,ytit=ytit3,yra=[0,0.48],ysty=1,$
        xra=[45.5,47.2],xsty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.1,0.1,0.32,0.39],/noerase;,xtickname=['44',' ','45',' ','46',' ','47']		
	frac=data[varlcuts].n_obsc / (data[varlcuts].n_obsc + data[varlcuts].n_unob)
	oplot,data[varlcuts].lcut,frac,linestyle=0,color=cgcolor('black'),thick=8



	xtit=textoidl('M_{switch}')
	plot,[0],[0],psym=0,$
        yra=[11.8,13.5],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.32,0.67,0.54,0.95],/noerase,xtickname=[' ',' ',' ',' ',' ',' '],$
        ytickname=[' ',' ',' ',' ',' ']					   
	oplot,data[varratio].ratio_lim,data[varratio].mhalo_obsc_med,linestyle=0,color=cgcolor('red'),thick=8
	oplot,data[varratio].ratio_lim,data[varratio].mhalo_unob_med,linestyle=0,color=cgcolor('blue'),thick=8
	oplot,data[varratio].ratio_lim,data[varratio].mhalo_obsc_mean,linestyle=2,color=cgcolor('red'),thick=8
	oplot,data[varratio].ratio_lim,data[varratio].mhalo_unob_mean,linestyle=2,color=cgcolor('blue'),thick=8
	
	plot,[0],[0],psym=0,$
        yra=[1.1,3.9],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.32,0.39,0.54,0.67],/noerase,xtickname=[' ',' ',' ',' ',' ',' '], $
        ytickname=[' ',' ',' ',' ',' ',' ',' ']		
	oplot,data[varratio].ratio_lim,data[varratio].b_obsc,linestyle=0,color=cgcolor('red'),thick=8
	oplot,data[varratio].ratio_lim,data[varratio].b_unob,linestyle=0,color=cgcolor('blue'),thick=8
	
	plot,[0],[0],psym=0,$
        xtit=xtit,yra=[0,0.48],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.32,0.1,0.54,0.39],/noerase,$
        ytickname=[' ',' ',' ',' ',' ',' '],xtickname=[' ','0.2','0.4','0.6','0.8',' ']					
	frac=data[varratio].n_obsc / (data[varratio].n_obsc + data[varratio].n_unob)
	oplot,data[varratio].ratio_lim,frac,linestyle=0,color=cgcolor('black'),thick=8



	xtit=textoidl('t_{obsc}/t_{tot}')
	plot,[0],[0],psym=0,$
        yra=[11.8,13.5],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.54,0.67,0.76,0.95],/noerase,xtickname=[' ',' ',' ',' ',' ',' '],	$
        ytickname=[' ',' ',' ',' ',' ']					   	   
	oplot,data[varobsclife].obsc_life_frac,data[varobsclife].mhalo_obsc_med,linestyle=0,color=cgcolor('red'),thick=8
	oplot,data[varobsclife].obsc_life_frac,data[varobsclife].mhalo_unob_med,linestyle=0,color=cgcolor('blue'),thick=8
	oplot,data[varobsclife].obsc_life_frac,data[varobsclife].mhalo_obsc_mean,linestyle=2,color=cgcolor('red'),thick=8
	oplot,data[varobsclife].obsc_life_frac,data[varobsclife].mhalo_unob_mean,linestyle=2,color=cgcolor('blue'),thick=8
	
	plot,[0],[0],psym=0,$
        yra=[1.1,3.9],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.54,0.39,0.76,0.67],/noerase,xtickname=[' ',' ',' ',' ',' ',' '],	$
        ytickname=[' ',' ',' ',' ',' ',' ',' ']			
	oplot,data[varobsclife].obsc_life_frac,data[varobsclife].b_obsc,linestyle=0,color=cgcolor('red'),thick=8
	oplot,data[varobsclife].obsc_life_frac,data[varobsclife].b_unob,linestyle=0,color=cgcolor('blue'),thick=8
	
	plot,[0],[0],psym=0,$
        xtit=xtit,yra=[0,0.48],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.54,0.1,0.76,0.39],/noerase,	$
        ytickname=[' ',' ',' ',' ',' ',' '],xtickname=[' ','0.2','0.4','0.6','0.8',' ']						
	frac=data[varobsclife].n_obsc / (data[varobsclife].n_obsc + data[varobsclife].n_unob)
	oplot,data[varobsclife].obsc_life_frac,frac,linestyle=0,color=cgcolor('black'),thick=8



	xtit=textoidl('Mean F_{edd}')
	plot,[0],[0],psym=0,$
        yra=[11.8,13.5],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.76,0.67,0.98,0.95],/noerase,xtickname=[' ',' ',' ',' ',' ',' '],	$
        ytickname=[' ',' ',' ',' ',' ']				   
	oplot,data[varfedd].mean_fedd,data[varfedd].mhalo_obsc_med,linestyle=0,color=cgcolor('red'),thick=5
	oplot,data[varfedd].mean_fedd,data[varfedd].mhalo_unob_med,linestyle=0,color=cgcolor('blue'),thick=5
	oplot,data[varfedd].mean_fedd,data[varfedd].mhalo_obsc_mean,linestyle=2,color=cgcolor('red'),thick=5
	oplot,data[varfedd].mean_fedd,data[varfedd].mhalo_unob_mean,linestyle=2,color=cgcolor('blue'),thick=5
	
	plot,[0],[0],psym=0,$
        yra=[1.1,3.9],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.76,0.39,0.98,0.67],/noerase,xtickname=[' ',' ',' ',' ',' ',' '],	$
        ytickname=[' ',' ',' ',' ',' ',' ',' ']				
	oplot,data[varfedd].mean_fedd,data[varfedd].b_obsc,linestyle=0,color=cgcolor('red'),thick=8
	oplot,data[varfedd].mean_fedd,data[varfedd].b_unob,linestyle=0,color=cgcolor('blue'),thick=8
	
	plot,[0],[0],psym=0,$
        xtit=xtit,yra=[0,0.48],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.76,0.1,0.98,0.39],/noerase,	$
        ytickname=[' ',' ',' ',' ',' ',' '],xtickname=[' ','0.2','0.4','0.6','0.8','1.0']					
	frac=data[varfedd].n_obsc / (data[varfedd].n_obsc + data[varfedd].n_unob)
	oplot,data[varfedd].mean_fedd,frac,linestyle=0,color=cgcolor('black'),thick=8
	PS_End



	PS_start,filename='params_v_nspurts.eps',xsize=19,ysize=5.5,/encapsul
 	circsym,/fill
	xtit=textoidl('Log(L_{cut} [ergs/s])')
	ytit=textoidl('Max N_{qso phases}')
	plot,[0],[0],psym=0,$
        xtit=xtit,ytit=ytit,yra=[0.5,4.5],ysty=1,$
        xra=[43.8,47.2],xsty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.1,0.22,0.32,0.95],/noerase,xtickname=['44',' ','45',' ','46',' ','47']		
	oplot,data[varlcuts].lcut,data[varlcuts].max_nspurts,linestyle=0,color=cgcolor('black'),thick=8

	xtit=textoidl('M_{bh}/M_{bh,exp} limit')
	plot,[0],[0],psym=0,$
        xtit=xtit,yra=[0.5,4.5],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.32,0.22,0.54,0.95],/noerase,$
        ytickname=[' ',' ',' ',' ',' ',' '],xtickname=[' ','0.2','0.4','0.6','0.8',' ']					
	oplot,data[varratio].ratio_lim,data[varratio].max_nspurts,linestyle=0,color=cgcolor('black'),thick=8

	xtit=textoidl('t_{obsc}/t_{tot}')
	plot,[0],[0],psym=0,$
        xtit=xtit,yra=[0.5,4.5],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.54,0.22,0.76,0.95],/noerase,	$
        ytickname=[' ',' ',' ',' ',' ',' '],xtickname=[' ','0.2','0.4','0.6','0.8',' ']						
	oplot,data[varobsclife].obsc_life_frac,data[varobsclife].max_nspurts,linestyle=0,color=cgcolor('black'),thick=8

	xtit=textoidl('Mean F_{edd}')
	plot,[0],[0],psym=0,$
        xtit=xtit,yra=[0.5,4.5],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.76,0.22,0.98,0.95],/noerase,	$
        ytickname=[' ',' ',' ',' ',' ',' '],xtickname=[' ','0.2','0.4','0.6','0.8','1.0']					
	oplot,data[varfedd].mean_fedd,data[varfedd].max_nspurts,linestyle=0,color=cgcolor('black'),thick=8
	PS_End
	


	PS_start,filename='params_v_redshift.eps',xsize=19,ysize=5.5,/encapsul
 	circsym,/fill
	xtit=textoidl('Log(L_{cut} [ergs/s])')
	ytit=textoidl('redshift')
	plot,[0],[0],psym=0,$
        xtit=xtit,ytit=ytit,yra=[0.7,1.5],ysty=1,$
        xra=[43.8,47.2],xsty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.1,0.22,0.32,0.95],/noerase,xtickname=['44',' ','45',' ','46',' ','47']		
	oplot,data[varlcuts].lcut,data[varlcuts].z_obsc_med,linestyle=0,color=cgcolor('red'),thick=8
	oplot,data[varlcuts].lcut,data[varlcuts].z_obsc_mean,linestyle=2,color=cgcolor('red'),thick=8
	oplot,data[varlcuts].lcut,data[varlcuts].z_unob_med,linestyle=0,color=cgcolor('blue'),thick=8
	oplot,data[varlcuts].lcut,data[varlcuts].z_unob_mean,linestyle=2,color=cgcolor('blue'),thick=8
	legend,['Mean','Median'],psym=[6,8],box=0,/top,/left,charsize=2.,charthick=2.

	xtit=textoidl('M_{bh}/M_{bh,exp} limit')
	plot,[0],[0],psym=0,$
        xtit=xtit,yra=[0.7,1.5],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.32,0.22,0.54,0.95],/noerase,$
        ytickname=[' ',' ',' ',' ',' ',' '],xtickname=[' ','0.2','0.4','0.6','0.8',' ']					
	oplot,data[varratio].ratio_lim,data[varratio].z_obsc_med,linestyle=0,color=cgcolor('red'),thick=8
	oplot,data[varratio].ratio_lim,data[varratio].z_obsc_mean,linestyle=2,color=cgcolor('red'),thick=8
	oplot,data[varratio].ratio_lim,data[varratio].z_unob_med,linestyle=0,color=cgcolor('blue'),thick=8
	oplot,data[varratio].ratio_lim,data[varratio].z_unob_mean,linestyle=2,color=cgcolor('blue'),thick=8

	xtit=textoidl('t_{obsc}/t_{tot}')
	plot,[0],[0],psym=0,$
        xtit=xtit,yra=[0.7,1.5],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.54,0.22,0.76,0.95],/noerase,	$
        ytickname=[' ',' ',' ',' ',' ',' '],xtickname=[' ','0.2','0.4','0.6','0.8',' ']						
	oplot,data[varobsclife].obsc_life_frac,data[varobsclife].z_obsc_med,linestyle=0,color=cgcolor('red'),thick=8
	oplot,data[varobsclife].obsc_life_frac,data[varobsclife].z_obsc_mean,linestyle=2,color=cgcolor('red'),thick=8
	oplot,data[varobsclife].obsc_life_frac,data[varobsclife].z_unob_med,linestyle=0,color=cgcolor('blue'),thick=8
	oplot,data[varobsclife].obsc_life_frac,data[varobsclife].z_unob_mean,linestyle=2,color=cgcolor('blue'),thick=8

	xtit=textoidl('Mean F_{edd}')
	plot,[0],[0],psym=0,$
        xtit=xtit,yra=[0.7,1.5],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.76,0.22,0.98,0.95],/noerase,	$
        ytickname=[' ',' ',' ',' ',' ',' '],xtickname=[' ','0.2','0.4','0.6','0.8','1.0']					
	oplot,data[varfedd].mean_fedd,data[varfedd].z_obsc_med,linestyle=0,color=cgcolor('red'),thick=8
	oplot,data[varfedd].mean_fedd,data[varfedd].z_obsc_mean,linestyle=2,color=cgcolor('red'),thick=8
	oplot,data[varfedd].mean_fedd,data[varfedd].z_unob_med,linestyle=0,color=cgcolor('blue'),thick=8
	oplot,data[varfedd].mean_fedd,data[varfedd].z_unob_mean,linestyle=2,color=cgcolor('blue'),thick=8
	PS_End
	

	print,' '
	print,' '


	match=where(data.b_obsc GT 2.7 AND data.b_obsc LT 3.5 $
				AND data.b_unob GT 1.6 AND data.b_unob LT 1.8 $
				AND data.z_obsc_med GT 0.95 AND data.z_obsc_med LT 1.05 $
				AND data.z_unob_med GT 0.95 AND data.z_unob_med LT 1.05,cnt)
;				AND data.lcut GT 45.8 AND data.lcut LT 46.2, cnt)
	frac=(data[match].n_obsc)/(data[match].n_obsc+data[match].n_unob)
		print,'The parameter combinations with'
		print,'    2.5 < b_obsc < 3.5         '
		print,'    1.55 < b_unob < 1.85         '
		print,'    0.9 < z < 1.1            '
		print,'          are                  '
		print,' (there are ' +strtrim(n_elements(match),2) + ')'
		print,' ----------------------------- '
	IF (cnt NE 0) THEN BEGIN
		FOR i=0L,n_elements(match)-1 DO BEGIN
			print,'L_cut:           ' + strtrim(data[match[i]].lcut,2)
			print,'M_bh/M_halo lim: ' + strtrim(data[match[i]].ratio_lim,2)
			print,'t_obsc/t_tot:    ' + strtrim(data[match[i]].obsc_life_frac,2)
			print,'Mean f_edd:      ' + strtrim(data[match[i]].mean_fedd,2)
			print,'   with   '
			print,'Med obsc M_halo: ' + strtrim(data[match[i]].mhalo_obsc_med,2)
			print,'Med unob M_halo: ' + strtrim(data[match[i]].mhalo_unob_med,2)
			print,'b_obsc:          ' + strtrim(data[match[i]].b_obsc,2)
			print,'b_unob:          ' + strtrim(data[match[i]].b_unob,2)
			print,'Mean/Med obsc z:      ' + strtrim(data[match[i]].z_obsc_mean,2) + ' / ' + strtrim(data[match[i]].z_obsc_med,2)
			print,'Mean/Med unob z:      ' + strtrim(data[match[i]].z_unob_mean,2) + ' / ' + strtrim(data[match[i]].z_unob_med,2)
			print,'f_obsc:          ' + strtrim(frac[i],2)
			print,'Max N_qsophases: ' + strtrim(data[match[i]].max_nspurts,2)
			print,'------------------------------------------------------------------'
		ENDFOR
	ENDIF ELSE BEGIN
		print,'No parameter combinations produce those results'
	ENDELSE
	print,' '

	xx=where(frac EQ max(frac))
	print,'The parameters that maximize f_obsc at ' + strtrim(max(frac),2) + ' are:'
	print,'L_cut:           ' + strtrim(data[match[xx]].lcut,2)
	print,'M_bh/M_halo lim: ' + strtrim(data[match[xx]].ratio_lim,2)
	print,'t_obsc/t_tot:    ' + strtrim(data[match[xx]].obsc_life_frac,2)
	print,'Mean f_edd:      ' + strtrim(data[match[xx]].mean_fedd,2)
	print,'   with   '
	print,'Med obsc M_halo: ' + strtrim(data[match[xx]].mhalo_obsc_med,2)
	print,'Med unob M_halo: ' + strtrim(data[match[xx]].mhalo_unob_med,2)
	print,'b_obsc:          ' + strtrim(data[match[xx]].b_obsc,2)
	print,'b_unob:          ' + strtrim(data[match[xx]].b_unob,2)
	print,'Mean/Med obsc z:      ' + strtrim(data[match[xx]].z_obsc_mean,2) + ' / ' + strtrim(data[match[xx]].z_obsc_med,2)
	print,'Mean/Med unob z:      ' + strtrim(data[match[xx]].z_unob_mean,2) + ' / ' + strtrim(data[match[xx]].z_unob_med,2)
	print,'Max N_qsophases: ' + strtrim(data[match[xx]].max_nspurts,2)
	print,' '
	print,' '
	
	PS_start,filename='params_match_observed.eps',xsize=22,ysize=6,/encapsul
 	circsym,/fill
	xtit=textoidl('Log(L_{cut} [ergs/s])')
	ytit=textoidl('N')
	plot,[0],[0],psym=0,$
        xtit=xtit,ytit=ytit,yra=[0.,25],ysty=1,$
        xra=[43.8,47.2],xsty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.7,ycharsize=1.7,$
        position=[0.06,0.22,0.24,0.95],/noerase,xtickname=['44',' ','45',' ','46',' ','47']		
	plothist,data[match].lcut,bin=0.1,thick=10,/half,/over,color=cgcolor('dark grey'),linestyle=2

	xtit=textoidl('M_{bh}/M_{bh,exp} limit')
	plot,[0],[0],psym=0,$
        xtit=xtit,yra=[0,25],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.24,0.22,0.42,0.95],/noerase,$
        ytickname=[' ',' ',' ',' ',' ',' '],xtickname=[' ','0.2','0.4','0.6','0.8',' ']					
	plothist,data[match].ratio_lim,bin=0.1,thick=10,/over,/half,color=cgcolor('dark grey'),linestyle=2

	xtit=textoidl('t_{obsc}/t_{tot}')
	plot,[0],[0],psym=0,$
        xtit=xtit,yra=[0.,25],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.42,0.22,0.60,0.95],/noerase,	$
        ytickname=[' ',' ',' ',' ',' ',' '],xtickname=[' ','0.2','0.4','0.6','0.8',' ']						
	plothist,data[match].obsc_life_frac,bin=0.1,thick=10,/over,/half,color=cgcolor('dark grey'),linestyle=2

	xtit=textoidl('Mean F_{edd}')
	plot,[0],[0],psym=0,$
        xtit=xtit,yra=[0.,25],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.60,0.22,0.78,0.95],/noerase,	$
        ytickname=[' ',' ',' ',' ',' ',' '],xtickname=[' ','0.2','0.4','0.6','0.8',' ']					
	plothist,data[match].mean_fedd,bin=0.1,thick=10,/over,/half,color=cgcolor('dark grey'),linestyle=2

	xtit=textoidl('f_{obsc}')
	plot,[0],[0],psym=0,$
        xtit=xtit,yra=[0.,25],ysty=1,$
        thick=5,xthick=8,ythick=8,$
        charthick=1.0,xcharsize=1.8,ycharsize=1.8,$
        position=[0.78,0.22,0.96,0.95],/noerase,	$
        ytickname=[' ',' ',' ',' ',' ',' '],xtickname=[' ','0.2','0.4','0.6','0.8','1.0']					
	plothist,frac,bin=0.1,thick=10,/over,/half,color=cgcolor('dark grey'),linestyle=2
	PS_End
	
	

	
stop
END

