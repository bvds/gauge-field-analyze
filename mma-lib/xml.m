xmlFilter::usage =
"Use expr=Import[\"filename.xml\"] to load XML file.
To explore the XML tree, use xmlFilter[expr[[2]],{\"tag1\", ...}, True]
To extract something, use xmlFilter[expr[[2]],{\"tag1\", ...}]";
xmlFilter[expr_String, {}, _] := ToExpression[expr];
xmlFilter[expr_XMLElement, {}, flag_:False] :=
  If[flag, expr[[1]], expr];
xmlFilter[XMLElement[a_, tags_, branches_], {a_, b___}, flag_:False] :=
  Select[Map[xmlFilter[#, {b}, flag]&, branches], (# =!= Null)&];
xmlFilter[XMLElement[_, tags_, branches_], {__}, _:False] := Null
