include (PKG)

PKG(ChargeTemplatePos
    DEPENDS
        Geometry
        EvtNavigator
        SimEvent
        BufferMemMgr

        ${SNIPER_RootWriter_LIBRARY}
)
