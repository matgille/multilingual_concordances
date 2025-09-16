<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:tei="http://www.tei-c.org/ns/1.0" xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:math="http://www.w3.org/2005/xpath-functions/math" exclude-result-prefixes="xs math"
    version="3.0">

    <xsl:output method="text"/>

    <xsl:template match="/">
        <xsl:result-document href="../data/prologue/la/prologue.txt">
            <xsl:apply-templates select="descendant::tei:div[@type = 'prologue']"/>
        </xsl:result-document>
    </xsl:template>


    <xsl:template match="tei:lb[@break = 'yes']">
        <xsl:text> </xsl:text>
    </xsl:template>

    <xsl:template match="tei:note | tei:teiHeader | tei:fw"/>

</xsl:stylesheet>
