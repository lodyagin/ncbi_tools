asntool -m ..\asn\general.asn -l asngen.h
asntool -m ..\asn\seqset.asn -l asnsset.h
asntool -m ..\asn\seqfeat.asn -l asnfeat.h -w100
asntool -m ..\asn\seqtable.asn -l asntable.h
asntool -m ..\asn\seqalign.asn -l asnalign.h
asntool -m ..\asn\seqloc.asn -l asnloc.h
asntool -m ..\asn\seqres.asn -l asnres.h
asntool -m ..\asn\medline.asn -l asnmedli.h
asntool -m ..\asn\biblio.asn -l asnbibli.h
asntool -m ..\asn\seq.asn -l asnseq.h
asntool -m ..\asn\pub.asn -l asnpub.h
asntool -m ..\asn\access.asn -l asnacces.h
asntool -m ..\asn\asnpub.all -l allpub.h
asntool -m ..\asn\asn.all -l all.h -w100
asntool -m ..\cdromlib\cdrom.asn -l cdrom.h
asntool -m ..\asn\submit.asn -l asnsubmt.h
asntool -m ..\asn\seqcode.asn -l asncode.h
asntool -m ..\asn\seqblock.asn -l asnblock.h
asntool -m ..\asn\featdef.asn -l asnfdef.h
asntool -m ..\network\entrez\client\netentr.asn -l asnneten.h
asntool -m ..\network\medarch\client\mla.asn -l asnmla.h
asntool -m ..\network\taxonomy\client\taxon.asn -l asntaxon.h
asntool -m ..\network\id0arch\client\id0.asn -l asnid0.h
asntool -m ..\network\spell\client\spell.asn -l spell.h
asntool -m ..\asn\objprt.asn -l asnprt.h
asntool -m ..\biostruc\mmdb1.asn -l mmdb1.h -w100
asntool -m ..\biostruc\mmdb2.asn -l mmdb2.h -w100
asntool -m ..\biostruc\mmdb3.asn -l mmdb3.h -w100
asntool -m ..\biostruc\cn3d\cn3d.asn -l cn3d.h -w100
rem asntool -m ..\network\blast2\client\blast18.asn -l asnbl18.h
asntool -m ..\network\blast3\client\blstspc.asn -l blstspc.h
asntool -m ..\network\suggest\client\suggest.asn -l sugmap.h
asntool -m ..\network\taxon1\common\stdalone\taxon1.asn -l taxon1.h
asntool -m ..\asn\ncbimime.asn -l asnmime.h
asntool -m ..\asn\pubmed.asn -l asnpubme.h
asntool -m ..\asn\medlars.asn -l asnmdrs.h
asntool -m ..\asn\proj.asn -l asnproj.h
