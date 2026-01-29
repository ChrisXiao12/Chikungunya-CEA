##****************************************************************************************
## FLOW DIAGRAM
##****************************************************************************************

if(isTRUE(create_flow_diagram)){
  g <- grViz(
    "
    digraph chikungunya_model {    
      graph [
        fontname = 'Latin Modern Math',
        layout = neato,
        splines = false,
        bgcolor = white
      ]
    
      node [
        shape = ellipse,
        fontname = 'Latin Modern Math',
        fontsize = 10,
        fontcolor = white,
        style = filled,
        fixedsize = true,
        width = 1.2,
        height = 0.7
      ]
    
      edge [
        fontname = 'Latin Modern Math',
        fontsize = 10,
        fontcolor = black,
        labeldistance = 0,
        labelangle = 0,
        labelfloat = false,
        arrowsize = 0.5
      ]
    
      S  [label='Susceptible,\nVaccine Naïve\n(S)', fillcolor='#E8F5E9', fontcolor=black, pos='0,2!']
      E  [label='Exposed\n(E)', fillcolor='#FBC02D', fontcolor=black, pos='2,2!']
      I  [label='Infectious\n(I)', fillcolor='#D32F2F', pos='4,2!']
      R  [label='Recovered\n(R)', fillcolor='#2E7D32', pos='6,2!']
    
      V  [label='Vaccinated\n(V)', fillcolor='#1B5E20', pos='0,0!']
      SV [label='Susceptible,\nVaccine Exposed\n(SV)', fillcolor='#E8F5E9', fontcolor=black, pos='2,0!']
      C  [label='Chronic Disease\n(C)', fillcolor='#B71C1C', pos='5,0!']
      D  [label='Death\n(D)', fillcolor='#6A1B9A', pos='3.5,-1.5!']
    
      SB1 [shape=point, width=0, height=0, style=invis, pos='-1.0,2!']
      SB2 [shape=point, width=0, height=0, style=invis, pos='-1.0,-1.5!']
      SB3 [shape=point, width=0, height=0, style=invis, pos='7.0,2!']
      SB4 [shape=point, width=0, height=0, style=invis, pos='7.0,-1.5!']
    
      S -> E [label='', arrowsize=0.5]
      betaslab [
        label = <<I>&beta;</I>>,
        shape = plaintext,
        fontsize = 10,
        fontcolor = black, 
        fontname = 'Latin Modern Math',
        width = 0,
        height = 0,
        pos = '0.8,2.1!'
      ];
      betassub [
        label = <<I>S</I>>, 
        shape = plaintext, 
        fontsize = 6,
        fontcolor = black, 
        fontname = 'Latin Modern Math', 
        width = 0, 
        height = 0, 
        pos = '0.84,2.05!'
      ];
      probs [
        label = <<I>(I/N)</I>>, 
        shape = plaintext, 
        fontsize = 10,
        fontcolor = black, 
        fontname = 'Latin Modern Math', 
        width = 0, 
        height = 0, 
        pos = '1.02,2.1!'
      ];
      S -> betaslab [style=invis];
      betaslab -> betassub [style=invis];
      betassub -> probs [style=invis];
      probs -> E [style=invis];

      S -> V [label = '', arrowsize = 0.5];
      epsilonlab [
        label = <<I>&epsilon;</I>>,
        shape = plaintext,
        fontsize = 10,
        fontcolor = black, 
        fontname = 'Latin Modern Math',
        width = 0,
        height = 0,
        pos = '-0.3,1.15!'
      ];
      epsilonsub [
        label = <<I>v</I>>, 
        shape = plaintext, 
        fontsize = 6,
        fontcolor = black, 
        fontname = 'Latin Modern Math', 
        width = 0, 
        height = 0, 
        pos = '-0.25,1.09!'
      ];
      psilab [
        label = <<I>&psi;</I>>,
        shape = plaintext,
        fontsize = 10,
        fontcolor = black, 
        fontname = 'Latin Modern Math',
        width = 0,
        height = 0,
        pos = '-0.18,1.15!'
      ];
      psisub [
        label = <<I>v</I>>, 
        shape = plaintext, 
        fontsize = 6,
        fontcolor = black, 
        fontname = 'Latin Modern Math', 
        width = 0, 
        height = 0, 
        pos = '-0.11,1.09!'
      ];
      S -> epsilonlab [style = invis];
      epsilonlab -> epsilonsub [style = invis];
      epsilonsub -> psilab [style = invis];
      psilab -> psisub [style = invis];
      psisub -> V [style = invis];

      S -> SV [label= '', arrowsize = 0.5];
      epsilonsvlab [
        label = <<I>(1-&epsilon;</I>>,
        shape = plaintext,
        fontsize = 10,
        fontcolor = black, 
        fontname = 'Latin Modern Math',
        width = 0,
        height = 0,
        pos = '0.40,1.15!'
      ];
      epsilonsvsub [
        label = <<I>v</I>>, 
        shape = plaintext, 
        fontsize = 6,
        fontcolor = black, 
        fontname = 'Latin Modern Math', 
        width = 0, 
        height = 0, 
        pos = '0.53,1.09!'
      ];
      psisvlab [
        label = <<I>)&psi;</I>>,
        shape = plaintext,
        fontsize = 10,
        fontcolor = black, 
        fontname = 'Latin Modern Math',
        width = 0,
        height = 0,
        pos = '0.62,1.15!'
      ];
      psisvsub [
        label = <<I>v</I>>, 
        shape = plaintext, 
        fontsize = 6,
        fontcolor = black, 
        fontname = 'Latin Modern Math', 
        width = 0, 
        height = 0, 
        pos = '0.71,1.09!'
      ];
      S -> epsilonsvlab [style = invis];
      epsilonsvlab -> epsilonsvsub [style = invis];
      epsilonsvsub -> psisvlab [style = invis];
      psisvlab -> psisvsub [style = invis];
      psisvsub -> SV [style = invis];

      ##V -> SV [label = <<I>&omega;<SUB>v</SUB></I> >]
      V -> SV [label = '', arrowsize = 0.5]
      omegalab [
        label = <<I>&omega;</I>>,
        shape = plaintext,
        fontsize = 10,
        fontcolor = black, 
        fontname = 'Latin Modern Math',
        width = 0,
        height = 0,
        pos = '0.95,0.08!'
      ];
      omegavsub [
        label = <<I>v</I>>, 
        shape = plaintext, 
        fontsize = 6,
        fontcolor = black, 
        fontname = 'Latin Modern Math', 
        width = 0, 
        height = 0, 
        pos = '1.04,0.04!'
      ];
      V -> omegalab [style = invis];
      omegalab -> omegavsub [style = invis];
      omegavsub -> SV [style = invis];

      V -> D [label = <<I>&mu;</I>       >]


      ##SV -> E [label = <<I>&beta;<SUB>SV</SUB>(I/N)</I> >]

      SV -> E [label='', arrowsize=0.5]
      betasvlab [
        label = <<I>&beta;</I>>,
        shape = plaintext,
        fontsize = 10,
        fontcolor = black, 
        fontname = 'Latin Modern Math',
        width = 0,
        height = 0,
        pos = '1.5,1.15!'
      ];
      betasvsub [
        label = <<I>SV</I>>, 
        shape = plaintext, 
        fontsize = 6,
        fontcolor = black, 
        fontname = 'Latin Modern Math', 
        width = 0, 
        height = 0, 
        pos = '1.56,1.09!'
      ];
      probsv [
        label = <<I>(I/N)</I>>, 
        shape = plaintext, 
        fontsize = 10,
        fontcolor = black, 
        fontname = 'Latin Modern Math', 
        width = 0, 
        height = 0, 
        pos = '1.78,1.16!'
      ];
      SV -> betasvlab [style=invis];
      betasvlab -> betasvsub [style=invis];
      betasvsub -> probsv [style=invis];
      probsv -> E [style=invis];

      SV -> D [label = <<I>&mu;</I>   >]
      E -> I [label = <<I>&lambda;</I>>]
      E:se -> D [label = <<I>&mu;</I> >]

      I -> R [label = '', arrowsize=0.5]
      gammarlab [
        label = <<I>&gamma;</I>>, 
        shape = plaintext, 
        fontsize = 10,
        fontcolor = black, 
        fontname = 'Latin Modern Math', 
        width = 0, 
        height = 0, 
        pos = '4.97,2.1!'
      ];
      gammarsub [
        label = <<I>R</I>>, 
        shape = plaintext, 
        fontsize = 6,
        fontcolor = black, 
        fontname = 'Latin Modern Math', 
        width = 0, 
        height = 0, 
        pos = '5.00,2.05!'
      ];
      I -> gammarlab [style=invis];
      gammarlab -> gammarsub [style=invis];
      gammarsub -> R [style=invis];

      I -> C [label= '', arrowsize = 0.5];
      L_IC [
        label = <<I>&gamma;</I>>,
        shape = plaintext,
        fontsize = 10,
        fontcolor = black, 
        fontname = 'Latin Modern Math',
        width = 0,
        height = 0,
        pos = '4.27,1.15!'
      ];
      R_IC [
        label = <<I>C</I>>, 
        shape = plaintext, 
        fontsize = 6,
        fontcolor = black, 
        fontname = 'Latin Modern Math', 
        width = 0, 
        height = 0, 
        pos = '4.3,1.09!'
      ];
      I -> L_IC [style = invis];
      L_IC -> R_IC [style = invis];
      R_IC -> C [style = invis];

      I -> D [label = <<I>&mu;+&delta;</I>>]
      C -> R [label = <<I>&kappa;</I> >]
      C -> D [label = <<I>&mu;</I>>]
    
      S:w -> SB1 [dir = none]
      SB1 -> SB2 [dir = none]
      SB2 -> D [headport = w, label = <<I>&mu;</I>>]
    
      R:e -> SB3 [dir = none]
      SB3 -> SB4 [dir = none]
      SB4 -> D [headport = e, label = <<I>&mu;</I> >]
    }
    "
    )
    
    ## Convert to SVG, then PNG
    svg <- export_svg(g)
    rsvg_png(
      svg = charToRaw(svg), 
      file = file.path(figures_dir, "chikungunya_model.png"), 
      width = 1600
    )
}

##*****************************************************************************************************
## END OF SCRIPT
##***************************************************************************************************** 
